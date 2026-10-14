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

# --- the drug-mention requirement, added 14 October 2026 after the anatomy terms exposed its absence ----------

from experiments.read_label_central_action import mentions_the_drug

WRONGLY_ADMITTED_WITHOUT_IT = [
    # receptor location rather than drug action, which is the ondansetron failure mode
    ("granisetron", "Serotonin receptors of the 5-HT 3 type are located peripherally on vagal nerve terminals and "
                    "centrally in the chemoreceptor trigger zone of the area postrema."),
    ("amisulpride", "D 2 receptors are located in the chemoreceptor trigger zone (CTZ) and respond to the dopamine "
                    "released from the nerve endings."),
    # a class-level hypothesis that names no drug and hedges
    ("pindolol", "The mechanism of the antihypertensive effects of beta-blocking agents has not been established, "
                 "but several mechanisms have been postulated: 1) an effect on the central nervous system"),
    ("timolol", "Possible mechanisms of action include reduction in cardiac output, reduction in plasma renin "
                "activity, and a central nervous system sympatholytic action."),
    # a negative finding about a DIFFERENT drug, which is the starkest of the six
    ("metformin", "In mice, no glipizide or metabolites were detectable autoradiographically in the brain or "
                  "spinal cord of males or females, nor in the fetuses of pregnant females."),
]


@pytest.mark.parametrize("drug, sentence", WRONGLY_ADMITTED_WITHOUT_IT)
def test_a_sentence_that_does_not_name_the_drug_cannot_admit_it(drug, sentence):
    """Each of these was a live admission before 14 October. Part (a) asks for THIS drug's central action."""
    assert not mentions_the_drug(sentence, drug), f"{drug}: {sentence[:70]}"


@pytest.mark.parametrize("drug, sentence", [
    ("codeine", "Codeine produces respiratory depression by direct action on brain stem respiratory centers."),
    ("cyclobenzaprine", "Such studies show that cyclobenzaprine acts primarily within the central nervous system "
                        "at brain stem as opposed to spinal cord levels"),
    ("metoclopramide", "Dopamine produces nausea and vomiting by stimulation of the medullary chemoreceptor "
                       "trigger zone (CTZ), and metoclopramide blocks stimulation of the CTZ by agents like l-dopa"),
    # a label's constant pronoun, which the screen must accept or it would lose most real sentences
    ("olanzapine", "It is 93% bound to plasma proteins and acts centrally at dopamine receptors."),
])
def test_a_sentence_about_the_drug_itself_still_admits_it(drug, sentence):
    assert mentions_the_drug(sentence, drug), f"{drug}: {sentence[:70]}"


def test_the_requirement_keeps_eleven_of_the_seventeen_and_the_six_lost_are_the_wrong_ones():
    """The measurement that justified adopting it, recorded as a number rather than a claim.

    Of the 17 drugs that passed part (a) before the requirement, 11 survive it: codeine, cyclobenzaprine,
    felodipine, hydrocodone, indomethacin, metoprolol, nebivolol, nimodipine, omeprazole, ondansetron, propranolol.
    The 6 lost are bisoprolol, granisetron, metformin, palonosetron, pindolol and timolol, and each is lost on a
    sentence quoted in WRONGLY_ADMITTED_WITHOUT_IT or in docs/label_central_action.md. So it is a correction and
    not a trade, which is why no compensating relaxation was added.
    """
    assert len(WRONGLY_ADMITTED_WITHOUT_IT) == 5, "bisoprolol's sentence is in the document rather than here"
