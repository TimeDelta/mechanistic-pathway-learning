"""Which label sentences name a drug's plasma carrier, and which are refused.

Every refusal here is a reading the first pass of experiments/scope_label_plasma_binder_coverage.py got wrong over the
130 cached labels, so each test quotes the sentence that produced the error rather than a constructed one.
"""
from experiments.scope_label_plasma_binder_coverage import binders_named, carrier_of_drug, quoted

STUDY_DRUGS = frozenset({"lidocaine", "prilocaine", "phenytoin", "fosphenytoin", "oxcarbazepine", "triazolam"})


def label(sentences: list[str], set_id: str = "abc12345", effective_time: str = "20250101") -> dict:
    return {"set_id": set_id, "effective_time": effective_time, "sentences": sentences}


def test_a_plain_statement_names_the_carrier():
    assert binders_named("Tiagabine is 96% bound to human plasma proteins, mainly to serum albumin and "
                         "α1-acid glycoprotein.", "tiagabine", STUDY_DRUGS) == {"albumin", "orosomucoid"}


def test_a_negated_binding_claim_names_no_carrier():
    assert binders_named("Oxcarbazepine and MHD do not bind to alpha-1-acid glycoprotein.",
                         "oxcarbazepine", STUDY_DRUGS) == set()


def test_a_negated_displacement_names_no_carrier():
    assert binders_named("Distribution Extremely high concentrations of triazolam do not displace bilirubin bound to "
                         "human serum albumin in vitro.", "triazolam", STUDY_DRUGS) == set()


def test_another_study_drug_as_the_binding_subject_names_no_carrier():
    assert binders_named("At concentrations produced by application of lidocaine and prilocaine cream, lidocaine is "
                         "approximately 70% bound to plasma proteins, primarily alpha-1-acid glycoprotein.",
                         "prilocaine", STUDY_DRUGS) == set()


def test_the_same_sentence_names_the_carrier_of_the_drug_it_is_about():
    assert binders_named("At concentrations produced by application of lidocaine and prilocaine cream, lidocaine is "
                         "approximately 70% bound to plasma proteins, primarily alpha-1-acid glycoprotein.",
                         "lidocaine", STUDY_DRUGS) == {"orosomucoid"}


def test_hypoalbuminemia_is_a_dosing_condition_rather_than_a_carrier():
    assert binders_named("Because the fraction of unbound phenytoin is increased in patients with renal or hepatic "
                         "disease, or in those with hypoalbuminemia, the monitoring of phenytoin serum levels should "
                         "be based on the unbound fraction.", "phenytoin", STUDY_DRUGS) == set()


def test_a_sentence_about_the_binding_of_other_drugs_names_no_carrier():
    assert binders_named("Donepezil hydrochloride at concentrations of 0.3 to 10 micrograms/mL did not affect the "
                         "binding of furosemide, digoxin and warfarin to human albumin.",
                         "donepezil", STUDY_DRUGS) == set()


def test_a_metabolite_sentence_carries_the_reading_without_naming_the_drug():
    """Oxcarbazepine's only surviving sentence is about MHD, so it names albumin and is marked as not naming the drug."""
    reading = carrier_of_drug({"drug": "oxcarbazepine", "labels": [
        label(["Approximately 40% of MHD is bound to serum proteins, predominantly to albumin.",
               "Oxcarbazepine and MHD do not bind to alpha-1-acid glycoprotein."])]}, STUDY_DRUGS)
    assert reading["carrier"] == "albumin"
    assert reading["names_the_drug"] is False


def test_a_drug_whose_every_sentence_is_refused_takes_no_carrier():
    reading = carrier_of_drug({"drug": "triazolam", "labels": [
        label(["Distribution Extremely high concentrations of triazolam do not displace bilirubin bound to human "
               "serum albumin in vitro."])]}, STUDY_DRUGS)
    assert reading["carrier"] == ""
    assert "sentence" not in reading


def test_the_quoted_sentence_states_the_binding_rather_than_the_dose():
    """Fosphenytoin's longest candidate is a dosing heading; the reading must rest on the binding sentence."""
    reading = carrier_of_drug({"drug": "fosphenytoin", "labels": [
        label(["2.7 Dosing in Patients with Renal or Hepatic Impairment or Hypoalbuminemia Because the fraction of "
               "unbound phenytoin (the active metabolite of fosphenytoin sodium injection) is increased in patients "
               "with renal or hepatic disease, the monitoring of serum levels should be based on albumin.",
               "Distribution Fosphenytoin is extensively bound (95% to 99%) to human plasma proteins, primarily "
               "albumin."])]}, STUDY_DRUGS)
    assert reading["carrier"] == "albumin"
    assert reading["sentence"].startswith("Distribution Fosphenytoin is extensively bound")


def test_the_quote_covering_most_of_the_carriers_is_preferred():
    """Donepezil's carrier is both, so the sentence naming both is the one a reader should be able to check."""
    reading = carrier_of_drug({"drug": "donepezil", "labels": [
        label(["Donepezil is approximately 96% bound to human plasma proteins, mainly to albumins (about 75%) and "
               "alpha 1 - acid glycoprotein (about 21%).",
               "Similarly, the binding of donepezil hydrochloride to human albumin was unchanged."])]}, STUDY_DRUGS)
    assert reading["carrier"] == "both"
    assert reading["sentence"].startswith("Donepezil is approximately 96% bound")


def test_a_long_run_together_table_is_quoted_around_its_carrier():
    sentence = ("Pharmacokinetic Parameters of Riluzole Absorption Bioavailability (oral) Approximately 60% Dose "
                "Proportionality Linear over a dose range of 25 mg to 100 mg every 12 hours Food effect AUC down 20% "
                "and Cmax down 45% high fat meal Distribution Plasma Protein Binding 96% (Mainly to albumin and "
                "lipoproteins) Elimination half-life 12 hours")
    excerpt = quoted(sentence, width=120)
    assert "albumin" in excerpt
    assert len(excerpt) <= 126


def test_a_short_sentence_is_quoted_whole():
    assert quoted("Entacapone binds mainly to serum albumin.") == "Entacapone binds mainly to serum albumin."
