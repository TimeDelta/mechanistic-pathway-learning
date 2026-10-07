"""Tests for the laboratory-abnormality labels: HPO logical definitions, ChEBI bridging and class-term exclusion."""
import pandas as pd

from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import (
    ChemicalDefinition,
    human_gem_metabolites_by_chebi,
    map_definition_to_metabolites,
    parse_hpo_chemical_definitions,
    read_chebi_relations,
)


def hpo_class(hpo_number: str, pato: str, chebi_number: str, uberon: str) -> str:
    return f'''    <owl:Class rdf:about="http://purl.obolibrary.org/obo/HP_{hpo_number}">
        <owl:equivalentClass>
            <owl:Restriction>
                <owl:someValuesFrom><owl:Class><owl:intersectionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="http://purl.obolibrary.org/obo/{pato}"/>
                    <rdf:Description rdf:about="http://purl.obolibrary.org/obo/CHEBI_{chebi_number}"/>
                    <owl:someValuesFrom rdf:resource="http://purl.obolibrary.org/obo/{uberon}"/>
                    <rdf:Description rdf:about="http://purl.obolibrary.org/obo/PATO_0000460"/>
                </owl:intersectionOf></owl:Class></owl:someValuesFrom>
            </owl:Restriction>
        </owl:equivalentClass>
    </owl:Class>
'''


OWL_TEXT = (hpo_class("0004923", "PATO_0000470", "28044", "UBERON_0000178")  # hyperphenylalaninemia
            + hpo_class("0000001", "PATO_0001997", "16526", "UBERON_0001088")  # decreased urinary carbon dioxide
            + hpo_class("0000002", "PATO_0000070", "26333", "UBERON_0000178")  # abnormal circulating prostaglandin
            + '    <owl:Class rdf:about="http://purl.obolibrary.org/obo/HP_0000003">\n        <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/HP_0000118"/>\n    </owl:Class>\n')

CHEBI_OBO = """[Term]
id: CHEBI:28044
name: phenylalanine
relationship: RO:0018034 CHEBI:32504 ! is conjugate acid of phenylalaninate

[Term]
id: CHEBI:32504
name: phenylalaninate

[Term]
id: CHEBI:17295
name: L-phenylalanine
is_a: CHEBI:28044 ! phenylalanine
relationship: RO:0018036 CHEBI:58095 ! is tautomer of L-phenylalanine zwitterion

[Term]
id: CHEBI:16998
name: D-phenylalanine
is_a: CHEBI:28044 ! phenylalanine

[Term]
id: CHEBI:58095
name: L-phenylalanine zwitterion

[Term]
id: CHEBI:16526
name: carbon dioxide

[Term]
id: CHEBI:26333
name: prostaglandin

[Term]
id: CHEBI:15551
name: prostaglandin E2
is_a: CHEBI:26333 ! prostaglandin

[Term]
id: CHEBI:15553
name: prostaglandin F2alpha
is_a: CHEBI:26333 ! prostaglandin

[Term]
id: CHEBI:15552
name: prostaglandin D2
is_a: CHEBI:26333 ! prostaglandin

[Term]
id: CHEBI:15554
name: prostaglandin H2
is_a: CHEBI:26333 ! prostaglandin
"""

HUMAN_GEM_METABOLITES = pd.DataFrame({
    "metsNoComp": ["MAM_PHE", "MAM_DPHE", "MAM_CO2", "MAM_PGE2", "MAM_PGF2A", "MAM_PGD2", "MAM_PGH2"],
    "metChEBIID": ["CHEBI:58095", "CHEBI:16998", "CHEBI:16526", "CHEBI:15551", "CHEBI:15553", "CHEBI:15552", "CHEBI:15554"],
})


def test_definitions_give_direction_chemical_and_fluid() -> None:
    definitions = {definition.hpo_id: definition for definition in parse_hpo_chemical_definitions(OWL_TEXT)}
    assert set(definitions) == {"HP:0004923", "HP:0000001", "HP:0000002"}  # a class without a chemical definition is skipped
    assert definitions["HP:0004923"] == ChemicalDefinition("HP:0004923", ("CHEBI:28044",), 1, "blood")
    assert definitions["HP:0000001"].direction == -1 and definitions["HP:0000001"].fluid == "urine"
    assert definitions["HP:0000002"].direction == 0


def test_bridging_reaches_the_l_form_and_leaves_out_d_forms_and_class_terms(tmp_path) -> None:
    obo_path = tmp_path / "chebi.obo"
    obo_path.write_text(CHEBI_OBO)
    names, neighbours, children = read_chebi_relations(obo_path)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(HUMAN_GEM_METABOLITES)
    definitions = {definition.hpo_id: definition for definition in parse_hpo_chemical_definitions(OWL_TEXT)}
    # phenylalanine: no Human-GEM entry for it or its conjugate base; the L- child, through its zwitterion, is the measured form
    assert map_definition_to_metabolites(definitions["HP:0004923"], names, neighbours, children, metabolites_by_chebi) == ({"MAM_PHE"}, "specific_form")
    assert map_definition_to_metabolites(definitions["HP:0000001"], names, neighbours, children, metabolites_by_chebi) == ({"MAM_CO2"}, "equivalent")
    # prostaglandin is a class: four specific forms exceed the limit of three
    assert map_definition_to_metabolites(definitions["HP:0000002"], names, neighbours, children, metabolites_by_chebi) == (set(), "class_term")
    excluded = frozenset({"MAM_CO2"})
    assert map_definition_to_metabolites(definitions["HP:0000001"], names, neighbours, children, metabolites_by_chebi, excluded) == (set(), "unmapped")


def test_potassium_ion_maps_to_human_gem_potassium_by_manual_entry(tmp_path) -> None:
    obo_path = tmp_path / "chebi.obo"
    obo_path.write_text("[Term]\nid: CHEBI:29103\nname: potassium(1+)\n\n[Term]\nid: CHEBI:26216\nname: potassium atom\n")
    names, neighbours, children = read_chebi_relations(obo_path)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(pd.DataFrame({"metsNoComp": ["MAM02200"], "metChEBIID": ["CHEBI:26216"]}))
    hypokalemia = ChemicalDefinition("HP:0002900", ("CHEBI:29103",), -1, "blood")
    assert map_definition_to_metabolites(hypokalemia, names, neighbours, children, metabolites_by_chebi) == ({"MAM02200"}, "manual")
