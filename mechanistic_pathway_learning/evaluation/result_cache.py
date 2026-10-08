"""On-disk cache of deterministic intermediate results, so that a long scoring job resumes after its container is
reclaimed instead of starting over (the baselines' per-fold fits, the aggregation's per-run entries and paired
bootstraps).

A result is stored under the SHA-256 of a key that names every input it depends on: array contents (array_digest),
settings and a digest of the source files that compute it (source_digest). A changed input or changed code gives a new
key, so a stale result is never reused; a matching key returns the stored result, which equals what recomputing it
would give for every deterministic computation (fixed bootstrap seeds, a seeded walk). Writes go to a temporary file
and are renamed into place, so a job killed mid-write leaves no truncated entry. Values are pickled, which keeps
tuples, integer dict keys and NaN exactly; the cache only ever holds what this code wrote itself.
"""
from __future__ import annotations

import hashlib
import json
import os
import pickle
from collections.abc import Callable
from pathlib import Path
from typing import TypeVar

import numpy as np

T = TypeVar("T")


def array_digest(*arrays) -> str:
    """SHA-256 over the dtype, shape and bytes of each array (None counts as its own value)."""
    digest = hashlib.sha256()
    for array in arrays:
        if array is None:
            digest.update(b"none;")
            continue
        values = np.ascontiguousarray(array)
        digest.update(f"{values.dtype.str}{values.shape};".encode())
        digest.update(values.tobytes())
    return digest.hexdigest()


def sparse_digest(matrix) -> str:
    """Digest of a scipy sparse matrix in CSR form."""
    csr = matrix.tocsr()
    return array_digest(np.array(csr.shape), csr.data, csr.indices, csr.indptr)


def source_digest(*paths: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(Path(path) for path in paths):
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def cached_result(cache_directory: Path | None, key: dict, compute: Callable[[], T]) -> T:
    """Return the stored result for key, or compute, store and return it. cache_directory None computes every time."""
    if cache_directory is None:
        return compute()
    name = hashlib.sha256(json.dumps(key, sort_keys=True, default=str).encode()).hexdigest()
    path = Path(cache_directory) / f"{name}.pkl"
    if path.exists():
        with path.open("rb") as handle:
            return pickle.load(handle)
    value = compute()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, path)
    return value
