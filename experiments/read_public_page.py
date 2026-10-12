"""Read a public web page from the shell, for the case where the fetch tool would wait for the owner's approval.

Why this exists. In the Claude app the fetch tool reads an address only when it came from a search result or from a
message of the owner, and asks the owner about any other address; the request is withdrawn after five minutes without
an answer. The owner's ruling of 11 October 2026 (CLAUDE.md) makes the shell the standing fallback for such pages.
This script is that fallback with its limits written into code, so that a research pass cannot forget them:

  1. the site's robots file is read first, and a page it disallows for automated readers is left unread;
  2. a page that answers with a check that the visitor is human is reported as refused and not saved;
  3. an address that carries a personal identifier (an e-mail address) is not requested;
  4. every page read is saved with its address, time and SHA-256, in a directory that holds nothing else.

It does not log in, solve challenges, or use caches and mirrors of a page that refused it.

Usage:
  python experiments/read_public_page.py https://example.org/page.html
  python experiments/read_public_page.py --output-directory data/raw/public_pages --print-text URL [URL ...]
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import urllib.robotparser
from pathlib import Path

USER_AGENT = "mechanistic-pathway-learning/0.1 (research reading of public pages; no login)"
# The robots rules are checked for this script's own name, for every reader ("*") and for the names under which the
# fetch tool and its provider's crawler announce themselves, because a site that turns those away has said it does
# not want an assistant reading the page.
ROBOTS_USER_AGENTS_CHECKED = ("mechanistic-pathway-learning", "*", "Claude-User", "ClaudeBot", "anthropic-ai")
REQUEST_TIMEOUT_SECONDS = 60
SECONDS_BETWEEN_REQUESTS_TO_ONE_HOST = 2.0
CHALLENGE_STATUS_CODES = frozenset({401, 403, 407, 429, 503})
CHALLENGE_MARKERS = (
    "verify you are a human", "verify you are human", "are you a robot", "captcha", "altcha-widget", "cf-chl",
    "challenge-platform", "enable javascript and cookies to continue", "unusual traffic", "access denied",
    # proof-of-work gates that answer with status 200 (a research pass of 11 October 2026 found two saved as pages)
    "anubis", "making sure you're not a bot", "making sure you&#39;re not a bot", "checking your browser", "checking if the site connection is secure", "proof-of-work", "proof of work",
)
DEFAULT_OUTPUT_DIRECTORY = Path("data/raw/public_pages")


def reason_address_is_not_requested(address: str) -> str | None:
    """Why this address must not be requested, or None when it may be."""
    parsed_address = urllib.parse.urlsplit(address)
    if parsed_address.scheme != "https":
        return "only https addresses are read"
    if not parsed_address.hostname or "." not in parsed_address.hostname:
        return "the address has no public host name"
    if parsed_address.username or parsed_address.password:
        return "the address carries credentials"
    decoded_address = urllib.parse.unquote(address)
    if re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", decoded_address):
        return "the address carries an e-mail address, and no personal identifier goes into a request"
    return None


def robots_file_allows(robots_file_text: str, address: str, user_agents: tuple[str, ...] = ROBOTS_USER_AGENTS_CHECKED) -> bool:
    """True when the robots rules allow the address for every one of the user agents checked."""
    parser = urllib.robotparser.RobotFileParser()
    parser.parse(robots_file_text.splitlines())
    return all(parser.can_fetch(user_agent, address) for user_agent in user_agents)


def looks_like_robots_file(response_text: str) -> bool:
    """A robots file is plain text; a site that answers the request for it with a web page has not served one."""
    leading_text = response_text.lstrip()[:400].lower()
    return not (leading_text.startswith("<!doctype") or leading_text.startswith("<html") or "<head" in leading_text)


def looks_like_challenge(status_code: int, response_text: str) -> bool:
    """True when the answer is a refusal or a check that the visitor is human rather than the page asked for."""
    if status_code in CHALLENGE_STATUS_CODES:
        return True
    leading_text = response_text[:20000].lower()
    return len(response_text) < 20000 and any(marker in leading_text for marker in CHALLENGE_MARKERS)


def text_of_page(page_source: str) -> str:
    """The readable text of an HTML page: scripts and styles removed, tags dropped, entities decoded."""
    without_scripts = re.sub(r"<(script|style|noscript)\b.*?</\1>", " ", page_source, flags=re.S | re.I)
    without_tags = re.sub(r"<[^>]+>", " ", without_scripts)
    return re.sub(r"\s+", " ", html.unescape(without_tags)).strip()


def request_bytes(address: str) -> tuple[int, str, bytes]:
    """Status code, content type and body of one GET. An HTTP error status is returned, not raised."""
    request = urllib.request.Request(address, headers={"User-Agent": USER_AGENT, "Accept": "text/html,application/xhtml+xml,text/plain,application/pdf;q=0.9,*/*;q=0.5"})
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            return response.status, response.headers.get("Content-Type", ""), response.read()
    except urllib.error.HTTPError as error:
        return error.code, error.headers.get("Content-Type", "") if error.headers else "", error.read() or b""


class PublicPageReader:
    """Reads pages one at a time, with one robots file per host and a pause between requests to the same host."""

    def __init__(self, output_directory: Path) -> None:
        self.output_directory = Path(output_directory)
        self.robots_file_by_host: dict[str, str | None] = {}
        self.last_request_time_by_host: dict[str, float] = {}

    def _wait_for_host(self, host: str) -> None:
        elapsed_seconds = time.monotonic() - self.last_request_time_by_host.get(host, 0.0)
        if elapsed_seconds < SECONDS_BETWEEN_REQUESTS_TO_ONE_HOST:
            time.sleep(SECONDS_BETWEEN_REQUESTS_TO_ONE_HOST - elapsed_seconds)
        self.last_request_time_by_host[host] = time.monotonic()

    def robots_file_of_host(self, host: str) -> str | None:
        """The host's robots rules as text; an empty string when it has none; None when it would not serve them."""
        if host not in self.robots_file_by_host:
            self._wait_for_host(host)
            status_code, _, body = request_bytes(f"https://{host}/robots.txt")
            robots_file_text = body.decode("utf-8", errors="replace")
            if status_code == 404:
                self.robots_file_by_host[host] = ""
            elif status_code == 200 and looks_like_robots_file(robots_file_text):
                self.robots_file_by_host[host] = robots_file_text
            else:
                self.robots_file_by_host[host] = None
        return self.robots_file_by_host[host]

    def read(self, address: str) -> dict:
        """One record per address: whether it was read, why not when it was not, and where the page was saved."""
        record: dict = {"address": address, "read": False, "requested_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        refusal_reason = reason_address_is_not_requested(address)
        if refusal_reason:
            return record | {"reason": refusal_reason}
        host = urllib.parse.urlsplit(address).hostname or ""
        try:
            robots_file_text = self.robots_file_of_host(host)
        except Exception as error:  # noqa: BLE001 - a network failure is reported, not raised
            return record | {"reason": f"the robots file could not be requested ({error})"}
        if robots_file_text is None:
            return record | {"reason": "the site did not serve a robots file to an automated reader, which is taken as a refusal"}
        if not robots_file_allows(robots_file_text, address):
            return record | {"reason": "the site's robots file disallows this address for automated readers"}
        self._wait_for_host(host)
        try:
            status_code, content_type, body = request_bytes(address)
        except Exception as error:  # noqa: BLE001
            return record | {"reason": f"the request failed ({error})"}
        is_text = "pdf" not in content_type.lower() and not body.startswith(b"%PDF")
        response_text = body.decode("utf-8", errors="replace") if is_text else ""
        record |= {"status_code": status_code, "content_type": content_type, "bytes": len(body)}
        if looks_like_challenge(status_code, response_text):
            return record | {"reason": f"the site answered with a refusal or a check that the visitor is human (HTTP {status_code})"}
        if status_code != 200:
            return record | {"reason": f"HTTP {status_code}"}
        digest = hashlib.sha256(body).hexdigest()
        self.output_directory.mkdir(parents=True, exist_ok=True)
        extension = ".pdf" if not is_text else (".html" if "html" in content_type.lower() else ".txt")
        saved_path = self.output_directory / f"{digest[:16]}{extension}"
        saved_path.write_bytes(body)
        record |= {"read": True, "sha256": digest, "saved_as": str(saved_path)}
        saved_path.with_suffix(saved_path.suffix + ".json").write_text(json.dumps(record, indent=1))
        if is_text:
            record["text"] = text_of_page(response_text) if extension == ".html" else response_text
        return record


def main() -> int:
    argument_parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    argument_parser.add_argument("addresses", nargs="+")
    argument_parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    argument_parser.add_argument("--print-text", action="store_true", help="print the readable text of each page read")
    argument_parser.add_argument("--text-characters", type=int, default=6000, help="how much text --print-text shows per page")
    arguments = argument_parser.parse_args()
    reader = PublicPageReader(arguments.output_directory)
    pages_not_read = 0
    for address in arguments.addresses:
        record = reader.read(address)
        page_text = record.pop("text", None)
        print(json.dumps(record))
        if not record["read"]:
            pages_not_read += 1
        elif arguments.print_text and page_text is not None:
            print(page_text[: arguments.text_characters])
    return 1 if pages_not_read else 0


if __name__ == "__main__":
    sys.exit(main())
