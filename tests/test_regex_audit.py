"""A negative-case audit of every classifying regular expression in the project.

Written 14 October 2026 on the owner's instruction ("please add a bunch of regex tests (for all regexes, not just
these ones) and review them critically"), after three classification passes in one day reported counts from
patterns nobody had checked the matches of. The common failure was never a missing positive case; it was a missing
NEGATIVE case - a phrase that matches and should not. So each test below pairs what a pattern must catch with the
near-misses that must not match it.

The patterns are grouped by what a wrong match costs:
  - a SIGN pattern decides whether a perturbation goes in at +1 or -1, so an error inverts the physiology;
  - an ADMISSION pattern decides whether a drug or a value enters the study at all;
  - a PARSING pattern pulls a number or an identifier out of structured text, where an error is usually loud;
  - a MECHANICAL pattern strips or splits text, where an error is quiet but bounded.
"""
import re

import pytest

# --- sign patterns -------------------------------------------------------------------------------------------

from experiments.scope_off_target_binding import NEGATIVE_ACTIONS, POSITIVE_ACTIONS, action_sign


@pytest.mark.parametrize("action, expected_sign", [
    ("Antagonist", -1), ("Inhibitor", -1), ("Channel blocker", -1), ("Negative allosteric modulator", -1),
    ("Inverse agonist", -1),  # handled before the two patterns, because "agonist" matches the positive one
    ("Agonist", 1), ("Partial agonist", 1), ("Activator", 1), ("Channel opener", 1),
    ("Positive allosteric modulator", 1), ("Potentiator", 1), ("Releasing agent", 1),
])
def test_the_action_sign_of_each_wording_gtopdb_uses(action, expected_sign):
    assert action_sign(action)[0] == expected_sign, action


@pytest.mark.parametrize("action", ["Inactivator", "Inactivates the channel", "Causes inactivation"])
def test_an_inactivating_action_is_not_read_as_positive(action):
    """"activat" sits inside "inactivat", which the audit found read these three as positive."""
    sign, _ = action_sign(action)
    assert sign == -1, f"{action!r} was classified {sign}"


def test_slowing_inactivation_is_positive_for_a_stated_physiological_reason():
    """The one live wording in GtoPdb 2026.3, on 23 rows. A channel whose inactivation is slowed stays open
    longer and carries more current, so positive is right - but the old pattern reached it through "activat"
    matching inside "inactivation", which is the wrong route to the right answer."""
    assert action_sign("Slows inactivation")[0] == 1


def test_a_wording_matching_both_patterns_is_left_unclassified_rather_than_guessed():
    assert action_sign("Antagonist and partial agonist")[0] is None


def test_the_partial_flag_is_independent_of_the_sign():
    assert action_sign("Partial agonist") == (1, True)
    assert action_sign("Agonist") == (1, False)


# --- admission patterns --------------------------------------------------------------------------------------

from experiments.build_drug_exposure_table import BRAIN_NOT_CEREBROSPINAL_FLUID, FLUID_OR_REGIONAL_VALUE


@pytest.mark.parametrize("column_name", [
    "Kp,uu,brain", "Kp,uu,bbb", "unbound brain", "Kp,uu brain", "Kp,uu, whole brain", "Kp,uu,wb",
    "Kpuu,bbb", "Kp,uu,ISF", "Kpuu,ecf",
])
def test_a_whole_brain_unbound_partition_column_is_recognised(column_name):
    assert BRAIN_NOT_CEREBROSPINAL_FLUID.search(column_name), column_name


@pytest.mark.parametrize("column_name", [
    "Kp,uu,CSF", "cerebrospinal fluid concentration", "Kp,uu,ROI", "Kp,uu,striatum", "lumbar CSF",
    "cisterna magna", "ventricle",
])
def test_a_fluid_or_regional_value_is_excluded(column_name):
    assert FLUID_OR_REGIONAL_VALUE.search(column_name), column_name


@pytest.mark.parametrize("column_name", [
    "regional Kp,uu", "Kp,uu,ROI", "Kp,uu,striatum", "Kp,uu,str", "Kp,uu hippocampus", "cerebellum Kp,uu",
    "thalamic Kp,uu", "frontal cortex Kp,uu",
])
def test_a_value_measured_in_one_part_of_the_brain_is_excluded(column_name):
    """What this screen must catch is a value for one part of the brain rather than across the whole of it,
    stated as a regional qualifier, as a region of interest, or by naming the region."""
    assert FLUID_OR_REGIONAL_VALUE.search(column_name), column_name


@pytest.mark.parametrize("not_a_regional_value", [
    "extracellular region", "brain region", "Kp,uu,brain", "unbound brain", "Kp,uu,bbb", "Kp,uu,wb",
])
def test_the_word_region_in_another_sense_no_longer_excludes_a_column(not_a_regional_value):
    """The bare term `region` used to match any text holding it, "extracellular region" - a Reactome compartment -
    and "brain region" among them. The two senses are now separated, so a whole-brain column survives."""
    assert not FLUID_OR_REGIONAL_VALUE.search(not_a_regional_value), not_a_regional_value


from experiments.scope_label_plasma_binder_coverage import (
    ALBUMIN_PATTERN,
    NEGATED_BINDING,
    OROSOMUCOID_PATTERN,
    STATES_BINDING,
    binders_named,
)


@pytest.mark.parametrize("label_text", [
    "alpha-1-acid glycoprotein", "alpha 1 acid glycoprotein", "orosomucoid", "AAG", "AGP",
    "\u03b1-1-acid glycoprotein",
])
def test_orosomucoid_is_recognised_under_the_names_a_label_uses(label_text):
    assert OROSOMUCOID_PATTERN.search(label_text.lower()), label_text


def test_the_orosomucoid_pattern_does_not_catch_an_unrelated_glycoprotein():
    """`(alpha|a)\\s?-?\\s?1?\\s?-?\\s?glycoprotein` is loose enough to be worth a negative case."""
    assert not OROSOMUCOID_PATTERN.search("p-glycoprotein")
    assert not OROSOMUCOID_PATTERN.search("glycoprotein iib/iiia")


def test_albumin_is_matched_as_a_word_start_so_prealbumin_is_not_albumin():
    assert ALBUMIN_PATTERN.search("serum albumin")
    assert not ALBUMIN_PATTERN.search("prealbumin"), "prealbumin is transthyretin, a different carrier"


def test_a_negated_binding_statement_is_detected_within_its_window():
    assert NEGATED_BINDING.search("the drug is not bound to plasma proteins")
    assert NEGATED_BINDING.search("does not displace warfarin")


def test_the_negation_window_ends_at_a_sentence_boundary_rather_than_running_on():
    """The window is 40 characters and stops at a full stop or semicolon, so a negation cannot reach the next
    sentence's binding statement. Without that the pattern would negate a neighbouring claim."""
    assert not NEGATED_BINDING.search("The drug is not metabolised. It is bound to albumin")


VALPROATE_ALBUMIN_SENTENCE = ("Liver disease is also associated with decreased albumin concentrations and larger "
                              "unbound fractions (2 to 2.6 fold increase) of valproate.")


def test_an_unbound_fraction_sentence_must_keep_counting_as_a_binding_statement():
    r"""Why STATES_BINDING's optional "un" is load-bearing and must not be narrowed.

    A refinement was proposed on 14 October 2026: exclude "fraction" when the "un" prefix is present, so that a
    fraction-unbound figure stops reading as a binding statement. It was tested against the data and rejected,
    because valproate's albumin carrier rests on exactly that construction. The sentence below is the one recorded
    in docs/label_plasma_binders.md and docs/label_plasma_binders_v4.md as valproate's albumin evidence; the
    proposed pattern `\bunbound\b(?!\s+fractions?)` does not match it, so the change would have dropped a carrier
    the screen was written to find.

    The looseness is therefore not a compromise but a requirement, and this test is what pins it.
    """
    assert STATES_BINDING.search(VALPROATE_ALBUMIN_SENTENCE.lower())
    narrowed = re.compile(r"\b(?:bind|binds|bound|binding)\b(?!\s+fraction)|\bunbound\b(?!\s+fractions?)")
    assert not narrowed.search(VALPROATE_ALBUMIN_SENTENCE.lower()), (
        "the narrowed form must be shown to fail, which is the reason it was not adopted")


OTHER_DRUGS_FOR_THE_GUARD = frozenset({"naloxone", "warfarin", "acetaminophen"})


@pytest.mark.parametrize("sentence, drug", [
    # the real sentence the subject guard was built for, in buprenorphine's label
    ("naloxone is approximately 45% protein bound, primarily to albumin", "buprenorphine"),
    # the same attribution error phrased with an unbound fraction, which the guard could not see until 14 October
    ("the unbound fraction of naloxone is determined by albumin concentration", "buprenorphine"),
    ("the unbound fraction of warfarin was unchanged; serum albumin was monitored", "celecoxib"),
    # an endogenous ligand competing for the carrier is not the drug either
    ("monitor when unbound bilirubin is elevated and serum albumin is low", "ceftriaxone"),
])
def test_an_unbound_figure_for_another_substance_gives_this_drug_no_carrier(sentence, drug):
    """The owner's standard, 14 October 2026: "if it COULD give a wrong carrier, it needs fixed anyways regardless
    of whether an existing example is present". None of these four appears in the cached labels; the first does and
    was already guarded. The fix is in the subject guard rather than in STATES_BINDING, so valproate survives."""
    assert binders_named(sentence, drug, OTHER_DRUGS_FOR_THE_GUARD) == set(), sentence


@pytest.mark.parametrize("sentence, drug, expected_carrier", [
    (VALPROATE_ALBUMIN_SENTENCE, "valproate", "albumin"),
    ("the drug is approximately 95% bound to human serum albumin", "phenytoin", "albumin"),
])
def test_the_guard_does_not_take_a_carrier_away_from_the_drug_the_sentence_is_about(sentence, drug, expected_carrier):
    assert binders_named(sentence, drug, OTHER_DRUGS_FOR_THE_GUARD) == {expected_carrier}, sentence


def test_a_protein_named_without_any_word_for_binding_states_no_carrier():
    """The complement of the test above: the screen exists because "Serum albumin less than 30 g/L", a dose rule,
    had given a radioligand an albumin carrier on 11 October 2026."""
    assert not STATES_BINDING.search("serum albumin less than 30 g/l")


# --- the central action patterns, which already have their own module ----------------------------------------

from experiments.read_label_central_action import (
    CENTRAL_ACTION, READS_AS_AN_ADVERSE_EVENT, READS_AS_NOT_CENTRAL, SENTENCE,
)


@pytest.mark.parametrize("adverse_event_text", [
    "Central Nervous System: Headache (8.2%), asthenia (0.6%), dizziness (3.1%)",
    "milder gastrointestinal and CNS side effects",
    "CNS adverse reactions were reported",
])
def test_an_adverse_event_table_or_phrase_is_refused(adverse_event_text):
    """Health Canada's diltiazem monograph carries the first of these as a table heading, which is why it matters."""
    assert READS_AS_AN_ADVERSE_EVENT.search(adverse_event_text), adverse_event_text


@pytest.mark.parametrize("negative_text", [
    "Only negligible amounts of the drug crossed the blood-brain barrier",
    "only minimal amounts were detected in the brain",
    "crosses the blood brain barrier poorly",
    "minimal distribution across the blood-brain barrier",
    "only minimal quantities have been detected in the cerebrospinal fluid",
])
def test_a_statement_that_the_drug_does_not_reach_the_brain_is_refused(negative_text):
    assert READS_AS_NOT_CENTRAL.search(negative_text), negative_text


def test_the_sentence_splitter_does_not_split_an_abbreviation_or_a_decimal():
    """The splitter needs a following capital, so "e.g. the" survives; a decimal has no following space."""
    assert len(SENTENCE.split("The value was 0.06 in the fluid.")) == 1
    assert len(SENTENCE.split("First sentence. Second sentence.")) == 2


def test_the_sentence_splitter_splits_on_a_semicolon_which_is_a_known_choice():
    """A semicolon split can cut a clause in half, which shortens a sentence and can drop the drug's name from it.

    Asserted so the behaviour is deliberate rather than discovered: a label clause like "X is an antagonist; it
    acts centrally" becomes two sentences, and the second loses the drug name to the pronoun rule.
    """
    assert len(SENTENCE.split("X is an antagonist; It acts centrally.")) == 2


# --- parsing patterns ----------------------------------------------------------------------------------------

from mechanistic_pathway_learning.graph.reactome_import import CHEBI_PATTERN, UNIPROT_PATTERN
from mechanistic_pathway_learning.graph.node_descriptors import SPECIES_PATTERN


def test_the_uniprot_pattern_takes_the_accession_and_not_the_surrounding_url():
    found = UNIPROT_PATTERN.findall('<rdf:li rdf:resource="http://identifiers.org/uniprot:P35367"/>')
    assert found == ["P35367"]


def test_the_chebi_pattern_takes_the_number_only():
    assert CHEBI_PATTERN.findall('chebiId=CHEBI:17996') == ["17996"]


def test_the_species_charge_pattern_reads_a_negative_charge():
    sbml = '<species id="M_MAM01442c" fbc:charge="-1" />'
    assert SPECIES_PATTERN.findall(sbml) == [("MAM01442", "-1")]


from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import MESH_DESCRIPTOR_PATTERN as DESCRIPTOR_ID


@pytest.mark.parametrize("identifier, is_valid", [
    ("D001921", True), ("D000067877", True), ("D01921", False), ("C001921", False), ("D001921x", False),
])
def test_the_mesh_descriptor_pattern_accepts_only_a_descriptor_identifier(identifier, is_valid):
    assert bool(DESCRIPTOR_ID.match(identifier)) == is_valid, identifier


from experiments.extract_brenda_enzyme_regulation import UNIPROT_ACCESSION_PATTERN


@pytest.mark.parametrize("text, expected", [
    ("P35367", ["P35367"]), ("Q9Y6L7", ["Q9Y6L7"]), ("A0A087X1C5", ["A0A087X1C5"]),
])
def test_the_brenda_accession_pattern_reads_the_accession_forms_uniprot_issues(text, expected):
    found = [match.group(1) for match in UNIPROT_ACCESSION_PATTERN.finditer(text)]
    assert found == expected, text


@pytest.mark.parametrize("not_an_accession", ["ATP", "NADPH", "GABA", "CHEBI"])
def test_the_brenda_accession_pattern_does_not_read_a_metabolite_abbreviation_as_an_accession(not_an_accession):
    assert not UNIPROT_ACCESSION_PATTERN.search(not_an_accession), not_an_accession


# --- mechanical patterns -------------------------------------------------------------------------------------

from mechanistic_pathway_learning.graph.tag_currency_metabolites import (
    CHARGE_SUFFIX_PATTERN as CHARGE_SUFFIX,
    COMPARTMENT_SUFFIX_PATTERN as COMPARTMENT_SUFFIX,
)


def test_the_compartment_suffix_is_stripped_in_both_spellings():
    assert COMPARTMENT_SUFFIX.sub("", "glucose[c]") == "glucose"
    assert COMPARTMENT_SUFFIX.sub("", "glucose_c") == "glucose"


def test_the_charge_suffix_strips_a_trailing_sign_and_leaves_an_internal_one():
    assert CHARGE_SUFFIX.sub("", "Na+") == "Na"
    assert CHARGE_SUFFIX.sub("", "Ca2+") == "Ca2"
    assert CHARGE_SUFFIX.sub("", "H2O") == "H2O"
