"""The guard that three classification passes on 14 October 2026 needed and did not have.

Each of those passes shipped a pattern without checking what it matched: a section comparison with no denominator,
a "mentions the CNS" screen where "medulla" caught the adrenal medulla and "neuron" caught sympathetic neurons, and
a widened action vocabulary where "medullary" caught the renal medullary thick ascending loop in celecoxib's label.
The first test below is the one that would have caught the last two.
"""
import re

import pytest

from experiments.read_label_central_action import (
    CENTRAL_ACTION,
    NOT_CENTRAL_DESPITE_THE_WORD,
    READS_AS_NOT_CENTRAL,
    vocabulary_does_not_match_peripheral_anatomy,
)


def test_the_vocabulary_matches_no_peripheral_anatomy():
    assert vocabulary_does_not_match_peripheral_anatomy(CENTRAL_ACTION) == []


def test_every_phrase_in_the_guard_list_is_actually_peripheral():
    """A guard list that drifted into holding central phrases would pass the test above for the wrong reason."""
    for phrase in NOT_CENTRAL_DESPITE_THE_WORD:
        assert not re.search(r"\b(?:brain|cerebral|cerebrospinal|central nervous)\b", phrase, re.I), phrase


def test_the_anatomy_terms_admit_the_sentence_they_were_added_for():
    metoclopramide_sentence = ("Dopamine produces nausea and vomiting by stimulation of the medullary chemoreceptor "
                               "trigger zone (CTZ), and metoclopramide blocks stimulation of the CTZ by agents like "
                               "l-dopa or apomorphine")
    assert CENTRAL_ACTION.search(metoclopramide_sentence)
    assert not READS_AS_NOT_CENTRAL.search(metoclopramide_sentence)


@pytest.mark.parametrize("penetration_sentence", [
    "steady-state predose concentrations of gabapentin in cerebrospinal fluid were approximately 20% of the corresponding plasma concentrations",
    "Theophylline passes freely across the placenta, into breast milk and into the cerebrospinal fluid (CSF)",
    "verapamil and norverapamil levels were noted in the cerebrospinal fluid with estimated partition coefficient of 0.06",
])
def test_a_penetration_statement_is_not_admitted_as_an_action(penetration_sentence):
    """Item 31 tightened part (a) from penetration to action, so a distribution figure must not admit a drug.

    These three sentences are what a wider vocabulary admitted when it was tried, which is why it was rejected.
    """
    assert not CENTRAL_ACTION.search(penetration_sentence), penetration_sentence
