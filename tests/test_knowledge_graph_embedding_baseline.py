"""The TransE baseline separates two chains that lead to different symptoms on a toy graph."""
import numpy as np

from mechanistic_pathway_learning.models.baselines.knowledge_graph_embedding_baseline import KnowledgeGraphEmbeddingBaseline


def test_transe_baseline_learns_chain_specific_symptoms() -> None:
    # two chains of 10 nodes with a shared hub (node 20) touching both; chain A genes induce symptom 0, chain B genes symptom 1
    edges = [(i, i + 1) for i in range(0, 9)] + [(i, i + 1) for i in range(10, 19)] + [(20, 0), (20, 10)]
    source = np.array([e[0] for e in edges])
    target = np.array([e[1] for e in edges])
    relation = np.array([0] * 18 + [1, 1])
    seeds = [np.array([i]) for i in range(20)]
    outcomes = np.zeros((20, 2))
    outcomes[:10, 0] = 1.0
    outcomes[10:, 1] = 1.0
    train = [i for i in range(20) if i % 3 != 0]
    test = [i for i in range(20) if i % 3 == 0]
    model = KnowledgeGraphEmbeddingBaseline(num_epochs=300, batch_size=64, embedding_dim=16, random_seed=0, device="cpu")
    model.fit(21, source, target, relation, 2, [seeds[i] for i in train], outcomes[train])
    predictions = model.predict([seeds[i] for i in test])
    assert predictions.shape == (len(test), 2)
    accuracy = np.mean([(predictions[k, 0] > predictions[k, 1]) == (outcomes[i, 0] == 1) for k, i in enumerate(test)])
    assert accuracy >= 0.7
    assert model.training_loss_trace[-1] < model.training_loss_trace[0]
    assert model.predict([np.array([], dtype=int)]).tolist() == [[0.0, 0.0]]  # a perturbation with no node scores zero
