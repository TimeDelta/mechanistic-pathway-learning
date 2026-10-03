"""Association types behind HPO annotations decide whether a monogenic pair can be grade A (review v0.4, finding 4)."""
from __future__ import annotations

from pathlib import Path

from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord, assign_evidence_grade
from mechanistic_pathway_learning.evidence.gene_disease_association_types import (
    association_type_for,
    load_omim_association_types,
    load_orphanet_association_types,
)


def test_lookup_tables_and_causality(tmp_path: Path) -> None:
    genes_to_disease = tmp_path / "genes_to_disease.txt"
    genes_to_disease.write_text("ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n"
                                "NCBIGene:1\tHMBS\tMENDELIAN\tOMIM:176000\tmim2gene\n"
                                "NCBIGene:2\tCOMT\tPOLYGENIC\tOMIM:181500\tmim2gene\n"
                                "NCBIGene:3\tGLT8D1\tUNKNOWN\tORPHA:803\torphadata\n")
    product6 = tmp_path / "en_product6.xml"
    product6.write_text("""<?xml version="1.0" encoding="UTF-8"?><JDBOR><DisorderList><Disorder id="1"><OrphaCode>803</OrphaCode>
<DisorderGeneAssociationList><DisorderGeneAssociation><Gene id="9"><Symbol>GLT8D1</Symbol></Gene>
<DisorderGeneAssociationType id="1"><Name lang="en">Major susceptibility factor in</Name></DisorderGeneAssociationType></DisorderGeneAssociation>
<DisorderGeneAssociation><Gene id="10"><Symbol>SOD1</Symbol></Gene>
<DisorderGeneAssociationType id="2"><Name lang="en">Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType></DisorderGeneAssociation>
</DisorderGeneAssociationList></Disorder></DisorderList></JDBOR>""")
    omim = load_omim_association_types(genes_to_disease)
    orphanet = load_orphanet_association_types(product6)
    assert association_type_for("HMBS", "OMIM:176000", omim, orphanet) == ("MENDELIAN", True)
    assert association_type_for("COMT", "OMIM:181500", omim, orphanet) == ("POLYGENIC", False)
    assert association_type_for("GLT8D1", "ORPHA:803", omim, orphanet) == ("Major susceptibility factor in", False)
    assert association_type_for("SOD1", "ORPHA:803", omim, orphanet) == ("Disease-causing germline mutation(s) in", True)
    assert association_type_for("SOD1", "ORPHA:999", omim, orphanet) == ("unknown", None)  # neither table types the pair


def test_grade_a_requires_a_causal_association_when_types_are_known() -> None:
    def record(known: bool, causal: int) -> EvidenceRecord:
        return EvidenceRecord("COMT", "psychosis", "induces", "monogenic", has_omim_clinical_synopsis=True, omim_entry_count=1,
                              association_type_known=known, causal_association_count=causal)

    assert assign_evidence_grade(record(known=True, causal=0)) == "C"  # polygenic or susceptibility only
    assert assign_evidence_grade(record(known=True, causal=1)) == "A"
    assert assign_evidence_grade(record(known=False, causal=0)) == "A"  # untyped pairs keep the provenance rule
