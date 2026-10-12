"""Tests: the per-drug exposure table (experiments/build_drug_exposure_table.py), the hand readings of protein binding
from drug labels (experiments/read_label_protein_binding.py) and the occupancy formula
(experiments/measure_target_occupancy.py). The readings are checked against the label sentence they quote, so a
reading whose numbers the sentence does not print is refused."""
import pytest

from experiments.build_drug_exposure_table import number_or_none, pubchem_compound_id, unbracketed_range
from experiments.measure_target_occupancy import occupancy_at
from experiments.read_label_protein_binding import numbers_are_printed, parsed_readings, verified_reading


def test_a_therapeutic_range_is_the_one_printed_outside_brackets() -> None:
    assert unbracketed_range("(5-15) 10-20") == (10.0, 20.0) and unbracketed_range("10-20 (5-15)") == (10.0, 20.0)
    assert unbracketed_range("0.05–0.2") == (0.05, 0.2) and unbracketed_range("appr. 0.5") == (0.5, 0.5)
    assert unbracketed_range("10-20; 30-40") is None and unbracketed_range("see text") is None and unbracketed_range("(5-15)") is None
    assert number_or_none("<0,5") == 0.5 and number_or_none("") is None and number_or_none(float("nan")) is None
    # the labels' compound identifiers: CID1 marks a stereo-merged compound and CID0 a stereo-specific one; an RxNorm concept has none
    assert pubchem_compound_id("CID100003121") == 3121 and pubchem_compound_id("CID000444") == 444 and pubchem_compound_id("RXCUI:10239") is None


def test_occupancy_is_one_half_at_the_dissociation_constant() -> None:
    assert occupancy_at(1.0, 6.0) == pytest.approx(0.5)  # 1 micromolar against a constant of 1 micromolar
    assert occupancy_at(0.01, 9.0) == pytest.approx(10.0 / 11.0)  # 10 nanomolar against 1 nanomolar
    assert occupancy_at(0.001, 6.0) == pytest.approx(1.0 / 1001.0)


def test_a_reading_is_kept_only_when_its_label_sentence_prints_its_numbers() -> None:
    label_file = {"labels": [
        {"set_id": "older", "effective_time": "20240101", "sentences": ["Distribution Examplamine is 90% -96% bound to serum albumin."]},
        {"set_id": "newer", "effective_time": "20260101", "sentences": ["Examplamine is 90%-96% bound to serum albumin.", "Bioavailability is approximately 60%."]}]}
    reading = {"drug": "examplamine", "percent_bound_low": 90.0, "percent_bound_high": 96.0, "qualifier": "", "fragment": "Examplamine is 90%-96% bound", "note": ""}
    verified = verified_reading(reading, label_file)
    assert verified["label_set_id"] == "newer" and verified["labels_with_the_sentence"] == 2  # white space does not matter; the newest label is quoted
    assert verified["fraction_unbound"] == pytest.approx(0.07) and verified["fraction_unbound_is"] == "a value" and not verified["bound_99_percent_or_more"]
    # a reading with numbers the sentence does not print is refused, and so is a fragment no label holds
    assert verified_reading({**reading, "percent_bound_high": 98.0}, label_file) is None
    assert verified_reading({**reading, "fragment": "Examplamine is 60% bound"}, label_file) is None
    # "greater than 99 percent bound" bounds the fraction unbound from above
    at_least = verified_reading({**reading, "percent_bound_low": 99.0, "percent_bound_high": 99.0, "qualifier": "at least", "fragment": "greater than 99% bound"},
                                {"labels": [{"set_id": "s", "effective_time": "2026", "sentences": ["It is greater than 99% bound to plasma proteins."]}]})
    assert at_least["fraction_unbound"] == pytest.approx(0.01) and at_least["fraction_unbound_is"] == "an upper limit" and at_least["bound_99_percent_or_more"]
    # a label that prints the free fraction supports the bound percentages derived from it
    assert numbers_are_printed({"percent_bound_low": 81.5, "percent_bound_high": 90.0}, "the free fraction increases from approximately 10% at 40 mcg/mL to 18.5% at 130 mcg/mL")
    assert not numbers_are_printed({"percent_bound_low": 80.0, "percent_bound_high": 90.0}, "the free fraction increases from approximately 10% at 40 mcg/mL to 18.5% at 130 mcg/mL")


def test_the_hand_readings_are_well_formed() -> None:
    readings = parsed_readings()
    assert len(readings) == len({reading["drug"] for reading in readings})  # one reading per drug
    assert all(0 <= reading["percent_bound_low"] <= reading["percent_bound_high"] <= 100 for reading in readings)
    assert {reading["qualifier"] for reading in readings} == {"", "at least", "at most"}
