"""Rate-limited, disk-cached HTTP GET for the NCBI and NLM web services used by the literature fetches.

Why: PubTator3, E-utilities and the MeSH RDF lookup are shared public services. NCBI asks for at most three
requests per second without an API key, and a fetch of several thousand identifiers must survive a dropped
connection without re-querying what it already has. Every response is written to the cache directory under
the SHA-256 of its URL (plus the POST body when there is one) before it is returned, so a rerun reads the
cache and makes no request at all; the cache file records the URL and the fetch time for the pin.
"""
from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

USER_AGENT = "mechanistic-pathway-learning/0.1 (literature evidence fetch)"
NCBI_MAX_REQUESTS_PER_SECOND = 3.0


class CachedRateLimitedClient:
    """GET (or POST with a form body) with a minimum interval between requests, retries with exponential
    backoff and a JSON cache keyed by request. Text responses are stored verbatim with their metadata."""

    def __init__(self, cache_directory: Path, max_requests_per_second: float = NCBI_MAX_REQUESTS_PER_SECOND, retries: int = 5, timeout_seconds: int = 60) -> None:
        self.cache_directory = Path(cache_directory)
        self.cache_directory.mkdir(parents=True, exist_ok=True)
        self.minimum_interval_seconds = 1.0 / max_requests_per_second
        self.retries = retries
        self.timeout_seconds = timeout_seconds
        self.last_request_time = 0.0
        self.request_count = 0
        self.cache_hit_count = 0

    @staticmethod
    def cache_key(url: str, post_body: bytes | None) -> str:
        digest = hashlib.sha256(url.encode("utf-8"))
        if post_body is not None:
            digest.update(b"\x00")
            digest.update(post_body)
        return digest.hexdigest()

    def cache_path(self, url: str, post_body: bytes | None = None) -> Path:
        return self.cache_directory / f"{self.cache_key(url, post_body)}.json"

    def fetch_text(self, url: str, post_body: bytes | None = None, accept: str = "application/json") -> str:
        """The response body as text, from the cache when present. Raises after the retries are exhausted."""
        path = self.cache_path(url, post_body)
        if path.exists():
            self.cache_hit_count += 1
            return json.loads(path.read_text(encoding="utf-8"))["body"]
        body_text = self._request_with_backoff(url, post_body, accept)
        payload = {"url": url, "post_body": None if post_body is None else post_body.decode("utf-8", errors="replace"), "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "body": body_text}
        temporary_path = path.with_suffix(".json.part")
        temporary_path.write_text(json.dumps(payload), encoding="utf-8")
        temporary_path.replace(path)
        return body_text

    def fetch_json(self, url: str, post_body: bytes | None = None):
        return json.loads(self.fetch_text(url, post_body))

    def _request_with_backoff(self, url: str, post_body: bytes | None, accept: str) -> str:
        last_error: Exception | None = None
        for attempt in range(self.retries):
            self._wait_for_rate_limit()
            request = urllib.request.Request(url, data=post_body, headers={"Accept": accept, "User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    self.request_count += 1
                    return response.read().decode("utf-8")
            except urllib.error.HTTPError as error:
                last_error = error
                if error.code in (400, 404) and attempt == 0:
                    # a malformed or unknown identifier is a stable answer, not a transient failure
                    raise
            except Exception as error:  # noqa: BLE001 - network retry
                last_error = error
            if attempt + 1 < self.retries:  # no backoff sleep after the final failed attempt
                time.sleep(min(60.0, 2.0 ** attempt))
        raise RuntimeError(f"request failed after {self.retries} attempts: {url}") from last_error

    def _wait_for_rate_limit(self) -> None:
        elapsed = time.monotonic() - self.last_request_time
        if elapsed < self.minimum_interval_seconds:
            time.sleep(self.minimum_interval_seconds - elapsed)
        self.last_request_time = time.monotonic()
