"""Tests for the on-disk result cache that lets the baselines and the aggregation resume after a reclaimed container."""
import numpy as np

from mechanistic_pathway_learning.evaluation.result_cache import array_digest, cached_result


def test_cached_result_computes_once_per_key_and_keeps_values_exactly(tmp_path) -> None:
    calls = []

    def compute() -> dict:
        calls.append(1)
        return {"point": float("nan"), "bins": {3: (0.1, 0.2)}}

    first = cached_result(tmp_path, {"what": "test", "inputs": array_digest(np.arange(4))}, compute)
    second = cached_result(tmp_path, {"what": "test", "inputs": array_digest(np.arange(4))}, compute)
    assert len(calls) == 1
    assert np.isnan(second["point"]) and second["bins"] == {3: (0.1, 0.2)} and first["bins"] == second["bins"]
    cached_result(tmp_path, {"what": "test", "inputs": array_digest(np.arange(5))}, compute)  # other inputs, other key
    assert len(calls) == 2
    assert not list(tmp_path.glob("*.tmp"))


def test_array_digest_tells_dtype_shape_and_none_apart() -> None:
    assert array_digest(np.zeros(4, dtype=np.float32)) != array_digest(np.zeros(4, dtype=np.float64))
    assert array_digest(np.zeros((2, 2))) != array_digest(np.zeros(4))
    assert array_digest(None, np.zeros(2)) != array_digest(np.zeros(2), None)


def test_without_a_directory_every_call_computes(tmp_path) -> None:
    calls = []
    for _ in range(2):
        cached_result(None, {"what": "test"}, lambda: calls.append(1))
    assert len(calls) == 2
