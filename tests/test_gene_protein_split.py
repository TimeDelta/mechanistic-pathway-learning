"""Tests: the gene and protein split (graph/gene_protein_split.py, docs/gene_protein_split.md) reroutes every relation
onto the protein nodes except the target end of transcription, gives genes of one shared UniProt entry one protein node
and a gene with several entries one protein node per group of entries the data tell apart (edge ends going to the
groups their source names) or one node for all of them, resolves renamed symbols, leaves genes without an entry as they
were, keeps
the gene brain expression block and expression columns on the genes and puts the protein block on the protein nodes,
makes the loader seed a drug's proteins and a knockout's gene with the leakage groups of the merged graph, widens seed
masking to the encodes partners and leaves encodes out of the rewiring."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from experiments.run_main_model import rewiring_fixed_relations
from mechanistic_pathway_learning.evaluation.experiment_data import GENE_TO_PROTEIN_FILE, load_experiment_data
from mechanistic_pathway_learning.evaluation.negative_controls import fast_degree_preserving_rewiring
from mechanistic_pathway_learning.graph.gene_protein_split import (
    ENCODES_RELATION,
    PROTEIN_NODE_TYPE,
    cell_class_weights_for_split,
    entries_named_by_interaction_rows,
    entries_of_gene_symbol,
    entry_groups_of_gene,
    place_descriptors_on_split,
    protein_node_assignment,
    resolve_entries_of_gene_symbols,
    resolved_protein_descriptors_on_merged,
    split_graph,
)
from mechanistic_pathway_learning.models.descriptor_treatments import seed_unit_mask
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

RELATIONS = ["catalyzed_by", "activates", "regulates_transcription_of", "member_of"]
NODE_IDS = ["GENE:TF", "GENE:B", "GENE:C", "GENE:D", "GENE:M", "R1", "COMPLEX:X"]
UNIPROT_ENTRIES = pd.DataFrame({"Entry": ["PTF", "PBC", "PD", "PM1", "PM2", "POTHER"],
                                "Gene Names (primary)": ["TF", "B; C", "D", "M", "M", "OTHER"]})
NAMED_ENTRIES = {("GENE:M", "GENE:D", "activates", "OmniPath", "source"): {"PM1"}}  # the M -> D rows name only PM1
DRUG_TARGET_ENTRIES = {"M": {"PM2"}}  # drug targets name PM2 for M
RECORDS_OF_ENTRY = {"PM1": {("edge end", "GENE:M", "GENE:D", "activates", "OmniPath", "source")}, "PM2": {("ChEMBL target", "CHEMBL1")}}  # the data tell M's entries apart


def toy_graph() -> tuple[pd.DataFrame, pd.DataFrame]:
    """TF regulates B; TF activates C and M; B and C (one shared entry) both activate D; M (two entries) activates D;
    D catalyses R1; D is in complex X."""
    nodes = pd.DataFrame({
        "node_id": NODE_IDS,
        "node_type": ["gene", "gene", "gene", "gene", "gene", "reaction", "protein_entity"],
        "gene_symbol": ["TF", "B", "C", "D", "M", None, None],
        "ensembl_gene_id": [None, "ENSG_B", None, "ENSG_D", "ENSG_M", None, None],
        "is_currency": [False] * 7, "degree": [0] * 7,
        "brain_median_tpm_max": [5.0, 2.0, 7.0, 1.0, 4.0, 3.0, None],
        "brain_expressed": [True, True, True, False, True, True, None]})
    edges = pd.DataFrame({
        "source_id": ["GENE:TF", "GENE:TF", "GENE:TF", "GENE:B", "GENE:C", "GENE:M", "GENE:D", "GENE:D"],
        "target_id": ["GENE:B", "GENE:C", "GENE:M", "GENE:D", "GENE:D", "GENE:D", "R1", "COMPLEX:X"],
        "relation_type": ["regulates_transcription_of", "activates", "activates", "activates", "activates", "activates", "catalyzed_by", "member_of"],
        "sign": [1.0] * 8, "evidence_source": ["CollecTRI", "OmniPath", "OmniPath", "OmniPath", "OmniPath", "OmniPath", "Human-GEM GPR", "Reactome entity member"]})
    return nodes, edges


def split_toy_graph(records_of_entry=RECORDS_OF_ENTRY):
    nodes, edges = toy_graph()
    assignment = protein_node_assignment(nodes, UNIPROT_ENTRIES, DRUG_TARGET_ENTRIES, records_of_entry)
    return nodes, edges, assignment, *split_graph(nodes, edges, RELATIONS, assignment, NAMED_ENTRIES)


def test_split_reroutes_relations_shares_one_entry_and_gives_several_entries_their_own_nodes() -> None:
    nodes, edges, assignment, split_nodes, split_edges, relations, summary = split_toy_graph()
    assert relations == RELATIONS + [ENCODES_RELATION]
    assert sorted(zip(assignment.gene_node_id, assignment.protein_node_id)) == [
        ("GENE:B", "PROTEIN:B"), ("GENE:C", "PROTEIN:B"), ("GENE:D", "PROTEIN:D"),
        ("GENE:M", "PROTEIN:M:PM1"), ("GENE:M", "PROTEIN:M:PM2"), ("GENE:TF", "PROTEIN:TF")]
    triples = set(zip(split_edges.source_id, split_edges.target_id, split_edges.relation_type))
    assert ("PROTEIN:TF", "GENE:B", "regulates_transcription_of") in triples  # the regulator's protein acts on the target's gene
    assert ("PROTEIN:TF", "PROTEIN:B", "activates") in triples
    assert ("PROTEIN:B", "PROTEIN:D", "activates") in triples  # B -> D and C -> D become one edge of the shared node
    assert ("PROTEIN:D", "R1", "catalyzed_by") in triples and ("PROTEIN:D", "COMPLEX:X", "member_of") in triples
    assert ("PROTEIN:M:PM1", "PROTEIN:D", "activates") in triples and ("PROTEIN:M:PM2", "PROTEIN:D", "activates") not in triples  # the entry the rows name
    assert {("PROTEIN:TF", "PROTEIN:M:PM1", "activates"), ("PROTEIN:TF", "PROTEIN:M:PM2", "activates")} <= triples  # no entry named: every entry
    assert summary["edge_ends_at_genes_with_several_proteins"] == {"named_subset": 1, "named_all": 0, "unresolved_all": 1}
    assert summary["genes_with_several_protein_nodes"] == 1 and summary["protein_nodes_of_genes_with_several"] == 2
    encodes = split_edges[split_edges.relation_type == ENCODES_RELATION]
    assert sorted(zip(encodes.source_id, encodes.target_id)) == sorted(zip(assignment.gene_node_id, assignment.protein_node_id))
    assert summary["edges_merged_after_rerouting"] == 1 and summary["shared_protein_nodes"] == 1
    assert set(split_edges[split_edges.source_id.str.startswith("GENE:")].relation_type) == {ENCODES_RELATION}  # a gene acts only through its proteins
    assert set(split_edges[split_edges.target_id.str.startswith("GENE:")].relation_type) == {"regulates_transcription_of"}
    by_node = split_nodes.set_index("node_id")
    assert by_node.degree["GENE:B"] == 2 and by_node.degree["GENE:C"] == 1 and by_node.degree["GENE:M"] == 2
    assert by_node.degree["PROTEIN:B"] == 4  # two encodes, the activation from TF, the activation of D
    protein = by_node.loc["PROTEIN:B"]
    assert protein.node_type == PROTEIN_NODE_TYPE and protein.display_name == "B/C" and protein.ensembl_gene_id == "ENSG_B" and protein.uniprot_entry == "PBC"
    assert by_node.loc["PROTEIN:M:PM2"].display_name == "M (PM2)"
    assert pd.isna(protein.brain_median_tpm_max) and not bool(protein.brain_expressed)  # expression stays on the genes (the user)
    assert by_node.loc["GENE:B"].brain_median_tpm_max == 2.0 and by_node.loc["R1"].brain_median_tpm_max == 3.0


def test_entries_the_data_do_not_tell_apart_share_one_protein_node() -> None:
    for records in (None, {"PM1": RECORDS_OF_ENTRY["PM1"]}):  # no record, or one entry named and the other never
        nodes, edges, assignment, split_nodes, split_edges, _, summary = split_toy_graph(records)
        assert sorted(zip(assignment.gene_node_id, assignment.protein_node_id))[-2:] == [("GENE:M", "PROTEIN:M"), ("GENE:TF", "PROTEIN:TF")]
        protein = assignment.set_index("protein_node_id").loc["PROTEIN:M"]
        assert protein.uniprot_entry == "PM1; PM2" and not protein.several_proteins and protein.drug_target
        triples = set(zip(split_edges.source_id, split_edges.target_id, split_edges.relation_type))
        assert {("PROTEIN:TF", "PROTEIN:M", "activates"), ("PROTEIN:M", "PROTEIN:D", "activates"), ("GENE:M", "PROTEIN:M", ENCODES_RELATION)} <= triples
        assert summary["genes_with_several_protein_nodes"] == 0 and summary["protein_nodes_of_several_entries"] == 1
        table = pd.DataFrame({"protein_rrr_1": [0.0] * 7, "protein_has_protein_descriptors": [0.0] * 7}, index=pd.Index(NODE_IDS, name="node_id"))
        per_entry = pd.DataFrame({"rrr_1": [0.5, 2.5, -1.0, 3.0, -2.0]}, index=pd.Index(["PTF", "PBC", "PD", "PM1", "PM2"], name="accession"))
        placed = place_descriptors_on_split(table, split_nodes, assignment, per_entry)
        assert placed.loc["PROTEIN:M", "protein_rrr_1"] == 0.5 and placed.loc["PROTEIN:M", "protein_has_protein_descriptors"] == 1.0  # the mean of its entries, the gene's row


def test_entry_groups_leave_out_entries_no_record_names_and_symbols_resolve_through_hgnc() -> None:
    records = {"A": {("edge end", 1)}, "B": {("edge end", 2), ("ChEMBL target", "T")}, "D": {("edge end", 2), ("ChEMBL target", "T")}}
    assert entry_groups_of_gene(["A", "B", "C", "D"], records) == ([["A"], ["B", "D"]], ["C"])  # B and D are named by the same records
    assert entry_groups_of_gene(["A", "C"], records) == ([["A", "C"]], [])  # one named group: no split
    assert entry_groups_of_gene([], records) == ([], [])
    uniprot = pd.DataFrame({"Entry": ["PNEW", "PACC", "PX"], "Gene Names (primary)": ["NEW", "ACC", "X"]})
    hgnc = pd.DataFrame({"symbol": ["NEW", "X", "Y"], "prev_symbol": ["OLD", "TWICE", "TWICE"]})
    entries, how = resolve_entries_of_gene_symbols(["NEW", "OLD", "PACC", "TWICE", "CHEBI:1"], uniprot, hgnc)
    assert entries == {"NEW": ["PNEW"], "OLD": ["PNEW"], "PACC": ["PACC"], "TWICE": [], "CHEBI:1": []}
    assert how == {"OLD": "previous HGNC symbol of NEW", "PACC": "reviewed accession"}
    nodes = pd.DataFrame({"node_id": ["GENE:NEW", "GENE:OLD", "GENE:CHEBI:1", "GENE:Z"], "node_type": ["gene"] * 4, "gene_symbol": ["NEW", "OLD", "CHEBI:1", "Z"],
                          "is_currency": [False] * 4, "degree": [0] * 4})
    edges = pd.DataFrame({"source_id": ["GENE:CHEBI:1", "GENE:OLD"], "target_id": ["GENE:Z", "GENE:CHEBI:1"], "relation_type": ["activates", "binds"],
                          "sign": [1.0, 1.0], "evidence_source": ["OmniPath", "OmniPath"]})
    uniprot = pd.concat([uniprot, pd.DataFrame({"Entry": ["PZ"], "Gene Names (primary)": ["Z"]})], ignore_index=True)
    assignment = protein_node_assignment(nodes, uniprot, hgnc_table=hgnc)
    assert sorted(zip(assignment.gene_node_id, assignment.protein_node_id)) == [("GENE:NEW", "PROTEIN:NEW"), ("GENE:OLD", "PROTEIN:NEW"), ("GENE:Z", "PROTEIN:Z")]
    _, split_edges, _, summary = split_graph(nodes, edges, RELATIONS, assignment)
    triples = set(zip(split_edges.source_id, split_edges.target_id, split_edges.relation_type))
    assert {("GENE:CHEBI:1", "PROTEIN:Z", "activates"), ("PROTEIN:NEW", "GENE:CHEBI:1", "binds")} <= triples  # the node without an entry keeps its edges
    assert summary["gene_nodes_without_a_protein_node"] == 1
    split_nodes = split_graph(nodes, edges, RELATIONS, assignment)[0]
    table = pd.DataFrame({"protein_rrr_1": [0.0, 0.0, 0.0, 4.0], "protein_has_protein_descriptors": [0.0, 0.0, 0.0, 1.0], "gene_brain_x": [1.0, 2.0, 3.0, 4.0]},
                         index=pd.Index(["GENE:NEW", "GENE:OLD", "GENE:CHEBI:1", "GENE:Z"], name="node_id"))
    per_entry = pd.DataFrame({"rrr_1": [7.0, 4.0]}, index=pd.Index(["PNEW", "PZ"], name="accession"))
    merged, filled = resolved_protein_descriptors_on_merged(table, assignment, place_descriptors_on_split(table, split_nodes, assignment, per_entry))
    assert filled == ["GENE:NEW", "GENE:OLD"] and merged.loc["GENE:OLD", "protein_rrr_1"] == 7.0 and merged.loc["GENE:Z", "protein_rrr_1"] == 4.0
    assert merged.loc["GENE:CHEBI:1", "protein_has_protein_descriptors"] == 0.0 and merged.gene_brain_x.tolist() == [1.0, 2.0, 3.0, 4.0]


def test_protein_block_sits_on_protein_nodes_one_row_per_entry_and_expression_stays_on_genes() -> None:
    nodes, _, assignment, split_nodes, _, _, _ = split_toy_graph()
    table = pd.DataFrame({"protein_rrr_1": [0.5, 0.0, 0.0, -1.0, 0.7, 0.0, 0.0], "protein_has_protein_descriptors": [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                          "gene_brain_region_x": [0.3, 1.2, -0.4, 0.0, 0.9, 0.0, 0.0], "reaction_brain_region_x": [0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0]},
                         index=pd.Index(NODE_IDS, name="node_id"))
    per_entry = pd.DataFrame({"rrr_1": [0.5, 2.5, -1.0, 3.0, -3.0]}, index=pd.Index(["PTF", "PBC", "PD", "PM1", "PM2"], name="accession"))
    placed = place_descriptors_on_split(table, split_nodes, assignment, per_entry)
    assert list(placed.index) == list(split_nodes.node_id)
    genes = ["GENE:TF", "GENE:B", "GENE:C", "GENE:D", "GENE:M"]
    assert (placed.loc[genes, ["protein_rrr_1", "protein_has_protein_descriptors"]] == 0).all().all()
    np.testing.assert_array_equal(placed.loc[genes, "gene_brain_region_x"], table.loc[genes, "gene_brain_region_x"])  # the genes keep their expression block
    proteins = split_nodes.node_id[split_nodes.node_type == PROTEIN_NODE_TYPE].tolist()
    assert (placed.loc[proteins, "gene_brain_region_x"] == 0).all()
    assert placed.loc["PROTEIN:B", "protein_rrr_1"] == 2.5 and placed.loc["PROTEIN:B", "protein_has_protein_descriptors"] == 1.0  # the shared entry's own row
    assert placed.loc["PROTEIN:M:PM1", "protein_rrr_1"] == 3.0 and placed.loc["PROTEIN:M:PM2", "protein_rrr_1"] == -3.0  # each entry its own row
    assert placed.loc["R1", "reaction_brain_region_x"] == 0.8 and placed.loc["PROTEIN:B", "reaction_brain_region_x"] == 0.0
    without_entries = place_descriptors_on_split(table, split_nodes, assignment)
    assert without_entries.loc["PROTEIN:M:PM1", "protein_rrr_1"] == 0.7 and without_entries.loc["PROTEIN:B", "protein_rrr_1"] == 0.0  # the gene rows' largest
    weights = pd.DataFrame({"neuron": [0.9, 0.2, 0.6, 0.1, 0.4, 0.5, 1.0]}, index=pd.Index(NODE_IDS, name="node_id"))
    split_weights = cell_class_weights_for_split(weights, split_nodes, assignment)
    assert split_weights.loc["GENE:B", "neuron"] == 0.2 and split_weights.loc["PROTEIN:B", "neuron"] == 0.6 and split_weights.loc["PROTEIN:TF", "neuron"] == 0.9
    assert split_weights.loc["PROTEIN:M:PM1", "neuron"] == 0.4 and split_weights.loc["PROTEIN:M:PM2", "neuron"] == 0.4
    assert not split_weights.isna().any().any()


def test_interaction_rows_name_the_entry_of_each_end_and_both_directions_of_binds() -> None:
    entries = entries_of_gene_symbol(UNIPROT_ENTRIES)
    assert entries["M"] == ["PM1", "PM2"] and entries["B"] == ["PBC"] and entries["C"] == ["PBC"] and entries["B; C"] == ["PBC"]
    rows = [{"source": "PM1", "target": "PD", "source_genesymbol": "M", "target_genesymbol": "D", "relation": "activates"},
            {"source": "COMPLEX:PM2_PTF", "target": "PD", "source_genesymbol": "M_TF", "target_genesymbol": "D", "relation": "binds"},
            # OmniPath sorts a complex's accessions but not its symbols: M takes PM1 here, not PTF in its position
            {"source": "PD", "target": "COMPLEX:PM1_PTF", "source_genesymbol": "D", "target_genesymbol": "TF_M", "relation": "inhibits"}]
    named = entries_named_by_interaction_rows(rows, "OmniPath", entries, lambda row: row["relation"])
    assert named == {("GENE:M", "GENE:D", "activates", "OmniPath", "source"): {"PM1"},
                     ("GENE:M", "GENE:D", "binds", "OmniPath", "source"): {"PM2"},
                     ("GENE:D", "GENE:M", "binds", "OmniPath", "target"): {"PM2"},
                     ("GENE:D", "GENE:M", "inhibits", "OmniPath", "target"): {"PM1"}}


def write_graph(directory: Path, nodes: pd.DataFrame, edges: pd.DataFrame, relations: list[str], assignment: pd.DataFrame | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    nodes.to_parquet(directory / "nodes.parquet", index=False)
    edges.to_parquet(directory / "edges.parquet", index=False)
    (directory / "relation_types.json").write_text(json.dumps(relations))
    if assignment is not None:
        assignment.to_parquet(directory / GENE_TO_PROTEIN_FILE, index=False)


def evidence_row(perturbation: str, kind: str, cluster: str, seeds: list) -> dict:
    return {"perturbation_id": perturbation, "perturbation_type": kind, "perturbation_label": perturbation, "group_id": perturbation,
            "disease_cluster_id": cluster, "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "A",
            "weight": 1.0, "label_frequency": None, "in_metabolic_layer": True, "evidence_date": None, "perturbation_nodes": json.dumps(seeds)}


def test_loader_seeds_a_drug_on_its_target_proteins_and_a_knockout_on_its_gene_with_the_merged_graph_groups(tmp_path: Path) -> None:
    nodes, edges, assignment, split_nodes, split_edges, relations, _ = split_toy_graph()
    write_graph(tmp_path / "merged", nodes, edges, RELATIONS)
    write_graph(tmp_path / "split", split_nodes, split_edges, relations, assignment)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row("knockout_D", "gene", "cluster:D", [["GENE:D", -1.0, 1.0]]),
                  evidence_row("knockout_M", "gene", "cluster:M", [["GENE:M", -1.0, 1.0]]),
                  evidence_row("drug_BC", "drug", "drug_BC", [["GENE:B", -1.0, 0.5], ["GENE:C", -1.0, 0.5], ["GENE:D", -1.0, 0.5]]),
                  evidence_row("drug_M", "drug", "drug_M", [["GENE:M", 1.0, 1.0]])]
                 ).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    merged = load_experiment_data(tmp_path / "merged", evidence_directory, group_by="disease_cluster_and_targets")
    split = load_experiment_data(tmp_path / "split", evidence_directory, group_by="disease_cluster_and_targets")
    seeds = {p: [split.node_ids[i] for i in s] for p, s in zip(split.perturbation_ids, split.perturbation_seeds)}
    assert seeds == {"drug_BC": ["PROTEIN:B", "PROTEIN:D"], "drug_M": ["PROTEIN:M:PM2"], "knockout_D": ["GENE:D"], "knockout_M": ["GENE:M"]}
    assert split.perturbation_magnitudes[split.perturbation_ids.index("drug_BC")].tolist() == [0.5, 0.5]
    assert split.group_ids == merged.group_ids  # each drug joins the knockout of its target on both graphs
    assert split.group_ids[split.perturbation_ids.index("drug_M")] == split.group_ids[split.perturbation_ids.index("knockout_M")]
    degree = dict(zip(split.perturbation_ids, split.perturbation_degrees))
    assert degree["knockout_M"] == 2.0  # the split degree of the seed itself (the user's decision): two encodes edges


def test_seed_masking_covers_the_encodes_partners_of_each_seed() -> None:
    perturbation_node_index = torch.tensor([[0, -1], [2, 4]])
    partners = torch.tensor([[0, 3, 2, 5], [3, 0, 5, 2]])  # 0 <-> 3 and 2 <-> 5
    assert seed_unit_mask(perturbation_node_index, 6).nonzero().tolist() == [[0, 0], [1, 2], [1, 4]]
    assert seed_unit_mask(perturbation_node_index, 6, partners).nonzero().tolist() == [[0, 0], [0, 3], [1, 2], [1, 4], [1, 5]]
    features = torch.randn(6, 4)
    without_partners = RelationalMessagePassingEncoder(6, 1, 8, 1, node_features=features, num_descriptor_columns=2, descriptor_treatment="seed_masked")
    encoder = RelationalMessagePassingEncoder(6, 1, 8, 1, node_features=features, num_descriptor_columns=2, descriptor_treatment="seed_masked",
                                              seed_mask_partner_index=partners)
    encoder.load_state_dict(without_partners.state_dict())  # the partner index is not saved, so checkpoints load either way
    sign_and_magnitude = torch.zeros(2, 2, 2)
    field = encoder.initial_node_state_field(perturbation_node_index, sign_and_magnitude, inject=False)
    torch.testing.assert_close(field[0, 3], encoder.base_state(without_descriptors=True)[3])  # the partner of seed 0 reads no descriptors
    torch.testing.assert_close(field[0, 1], encoder.base_state()[1])  # other nodes keep theirs
    field_without_partners = without_partners.initial_node_state_field(perturbation_node_index, sign_and_magnitude, inject=False)
    torch.testing.assert_close(field_without_partners[0, 3], without_partners.base_state()[3])


def test_rewiring_leaves_encodes_fixed_and_the_other_relations_as_without_it() -> None:
    generator = np.random.default_rng(3)
    activation_edges = generator.choice(40, size=(2, 120))
    activation_edges = activation_edges[:, activation_edges[0] != activation_edges[1]]
    encodes_edges = np.stack([np.arange(40, 60), np.arange(60, 80)])
    edges = np.concatenate([activation_edges, encodes_edges], axis=1)
    relations = np.array([0] * activation_edges.shape[1] + [1] * encodes_edges.shape[1])
    fixed = fast_degree_preserving_rewiring(edges, relations, num_swaps_per_edge=10, random_seed=7, fixed_relations=[1])
    np.testing.assert_array_equal(fixed[:, relations == 1], encodes_edges)
    free = fast_degree_preserving_rewiring(edges, relations, num_swaps_per_edge=10, random_seed=7)
    np.testing.assert_array_equal(fixed[:, relations == 0], free[:, relations == 0])  # encodes comes last, so the stream before it is the same
    assert not np.array_equal(free[:, relations == 1], encodes_edges)
    assert rewiring_fixed_relations(RELATIONS, False) == []
    assert rewiring_fixed_relations(RELATIONS + [ENCODES_RELATION], False) == [len(RELATIONS)]
    assert rewiring_fixed_relations(RELATIONS + [ENCODES_RELATION], True) == []
