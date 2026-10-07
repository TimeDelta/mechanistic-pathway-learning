import numpy as np

from experiments.score_by_module_coverage import largest_tie_group


def test_largest_tie_group_marks_the_rows_sharing_one_prediction():
    bias_row = [0.1, 0.2]
    scores = np.array([bias_row, [0.5, 0.6], [0.1 + 1e-8, 0.2], bias_row, [0.3, 0.3], [0.3, 0.3]])
    assert largest_tie_group(scores).tolist() == [True, False, True, True, False, False]


def test_largest_tie_group_is_empty_when_no_two_rows_tie():
    scores = np.array([[0.1, 0.2], [0.1, 0.3], [0.2, 0.2]])
    assert not largest_tie_group(scores).any()
    assert not largest_tie_group(scores[:1]).any()


def test_largest_tie_group_compares_every_symptom():
    # equal in the first symptom only, so not tied
    scores = np.array([[0.1, 0.2], [0.1, 0.4], [0.1, 0.6]])
    assert not largest_tie_group(scores).any()
