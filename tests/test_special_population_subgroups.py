"""Which special-population subgroups the carrier reading uses, and which leave on age.

The user, 10 October 2026: "I didn't say PEDIATRIC is irrelevant. Only NEWBORN through TODDLER. At the age of five,
per much everything should be diagnosable other than I.e. bipolar disorder". The rule is therefore a band of ages and
not a population, so these tests pin that an older child survives it, which the pinned database cannot show because
it holds no row above infant.
"""
import math

import pandas as pd

from experiments.measure_special_population_carrier_signal import (SUBGROUPS_BELOW_THE_DIAGNOSABLE_AGE, SUBGROUPS_TO_DECIDE,
                                                                   leave_one_drug_out_predictions, population_readings)

EXCLUDED = {name.lower() for name in SUBGROUPS_BELOW_THE_DIAGNOSABLE_AGE}


def rows_with_subgroups(subgroups: list[str]) -> pd.DataFrame:
    """Two carriers per subgroup, enough for a reading, with the ratio ordered so orosomucoid sits above albumin."""
    collected = []
    for index, subgroup in enumerate(subgroups):
        for carrier, ratio in (("albumin", 0.9), ("orosomucoid", 1.4)):
            collected.append({"carrier": carrier, "drug": f"{carrier}_drug_{index}", "subgroup": subgroup,
                              "population": "paediatric", "fraction_unbound_reference": 0.1,
                              "fraction_unbound_special": 0.1 * ratio, "log_ratio": math.log(ratio)})
    return pd.DataFrame(collected)


def kept(rows: pd.DataFrame) -> pd.DataFrame:
    return rows[~rows.subgroup.str.strip().str.lower().isin(EXCLUDED)]


def test_newborn_through_toddler_leaves():
    for subgroup in ("Neonate", "Newborn", "Preterm", "Infant", "Toddler"):
        assert subgroup.lower() in EXCLUDED, subgroup


def test_an_older_child_stays():
    for subgroup in ("Child", "School age", "Adolescent", "Teenager"):
        assert subgroup.lower() not in EXCLUDED, subgroup


def test_preschool_is_left_for_the_user_to_decide_rather_than_dropped():
    """The line sits between toddler and five, so a 3-to-5 band is neither excluded nor silently kept."""
    assert "preschool" not in EXCLUDED
    assert "preschool" in SUBGROUPS_TO_DECIDE


def test_a_school_age_row_survives_the_rule_while_an_infant_row_does_not():
    rows = rows_with_subgroups(["Neonate", "Infant", "School age"])
    assert set(kept(rows).subgroup) == {"School age"}


def test_a_reading_is_still_produced_from_the_subgroups_that_stay():
    rows = rows_with_subgroups(["Neonate", "School age"])
    readings = population_readings(kept(rows))
    assert len(readings) == 1
    assert int(readings.iloc[0].measurements) == 2


def test_dropping_every_subgroup_leaves_no_drug_classified():
    rows = rows_with_subgroups(["Neonate", "Infant"])
    assert len(leave_one_drug_out_predictions(kept(rows))) == 0
