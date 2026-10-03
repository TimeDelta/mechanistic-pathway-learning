"""Baseline B2, knowledge-graph embedding (design section 5.6).

The physiology graph and the training evidence are one typed knowledge graph: every graph edge is
a (head, relation, tail) triple, each symptom is a node, and every training observation adds a
(perturbation node, induces, symptom) triple. TransE [Bordes et al. 2013] or ComplEx embeddings are
trained with margin ranking against corrupted heads and tails sampled uniformly (the degree-matched
alternative of configs/baselines is left for the full pipeline). The score of a held-out perturbation
for a symptom is the plausibility of (perturbation node, induces, symptom); for a perturbation with
several seed nodes (a drug with a protein-family target) the scores of its nodes are averaged.

This is the "densely connected entities being highly ranked no matter the context" baseline of
assumption A9: hub nodes get strong embeddings because they appear in many triples, so its margin
over the popularity and degree baselines is what the model has to beat on top of hub structure.
Pure torch; CPU-sized (a few minutes for 25,000 nodes and 100,000 triples).
"""
from __future__ import annotations

import numpy as np
import torch
from torch import nn


class TransEScorer(nn.Module):
    def __init__(self, num_entities: int, num_relations: int, embedding_dim: int) -> None:
        super().__init__()
        self.entity = nn.Embedding(num_entities, embedding_dim)
        self.relation = nn.Embedding(num_relations, embedding_dim)
        bound = 6.0 / np.sqrt(embedding_dim)
        nn.init.uniform_(self.entity.weight, -bound, bound)
        nn.init.uniform_(self.relation.weight, -bound, bound)
        with torch.no_grad():
            self.relation.weight.div_(self.relation.weight.norm(dim=1, keepdim=True).clamp_min(1e-9))

    def forward(self, head: torch.Tensor, relation: torch.Tensor, tail: torch.Tensor) -> torch.Tensor:
        """Higher is more plausible: minus the L1 distance of head + relation from tail."""
        return -(self.entity(head) + self.relation(relation) - self.entity(tail)).abs().sum(dim=-1)

    def normalize_entities(self) -> None:
        with torch.no_grad():
            self.entity.weight.div_(self.entity.weight.norm(dim=1, keepdim=True).clamp_min(1.0))


class ComplExScorer(nn.Module):
    def __init__(self, num_entities: int, num_relations: int, embedding_dim: int) -> None:
        super().__init__()
        self.entity_real = nn.Embedding(num_entities, embedding_dim)
        self.entity_imaginary = nn.Embedding(num_entities, embedding_dim)
        self.relation_real = nn.Embedding(num_relations, embedding_dim)
        self.relation_imaginary = nn.Embedding(num_relations, embedding_dim)
        for table in (self.entity_real, self.entity_imaginary, self.relation_real, self.relation_imaginary):
            nn.init.normal_(table.weight, std=0.1)

    def forward(self, head: torch.Tensor, relation: torch.Tensor, tail: torch.Tensor) -> torch.Tensor:
        head_real, head_imaginary = self.entity_real(head), self.entity_imaginary(head)
        tail_real, tail_imaginary = self.entity_real(tail), self.entity_imaginary(tail)
        relation_real, relation_imaginary = self.relation_real(relation), self.relation_imaginary(relation)
        return (relation_real * head_real * tail_real + relation_real * head_imaginary * tail_imaginary
                + relation_imaginary * head_real * tail_imaginary - relation_imaginary * head_imaginary * tail_real).sum(dim=-1)

    def normalize_entities(self) -> None:
        return None


class KnowledgeGraphEmbeddingBaseline:
    def __init__(self, model_name: str = "TransE", embedding_dim: int = 64, num_epochs: int = 40, batch_size: int = 4096, learning_rate: float = 0.01,
                 margin: float = 1.0, negatives_per_positive: int = 4, random_seed: int = 0, device: str | None = None) -> None:
        self.model_name = model_name
        self.embedding_dim = embedding_dim
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.margin = margin
        self.negatives_per_positive = negatives_per_positive
        self.random_seed = random_seed
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.scorer: nn.Module | None = None
        self.num_graph_nodes = 0
        self.num_symptoms = 0
        self.induces_relation_index = 0
        self.training_loss_trace: list[float] = []

    def fit(self, num_graph_nodes: int, edge_source: np.ndarray, edge_target: np.ndarray, edge_relation: np.ndarray, num_relation_types: int,
            training_perturbation_seeds: list[np.ndarray], training_outcomes: np.ndarray) -> "KnowledgeGraphEmbeddingBaseline":
        """Graph triples plus (seed node, induces, symptom node) triples for every training positive."""
        torch.manual_seed(self.random_seed)
        generator = np.random.default_rng(self.random_seed)
        self.num_graph_nodes = num_graph_nodes
        self.num_symptoms = training_outcomes.shape[1]
        self.induces_relation_index = num_relation_types
        heads = [np.asarray(edge_source, dtype=np.int64)]
        relations = [np.asarray(edge_relation, dtype=np.int64)]
        tails = [np.asarray(edge_target, dtype=np.int64)]
        for seeds, outcomes in zip(training_perturbation_seeds, training_outcomes):
            for symptom_index in np.where(outcomes > 0)[0]:
                for seed in seeds:
                    heads.append(np.array([seed], dtype=np.int64))
                    relations.append(np.array([self.induces_relation_index], dtype=np.int64))
                    tails.append(np.array([num_graph_nodes + symptom_index], dtype=np.int64))
        head = torch.as_tensor(np.concatenate(heads), device=self.device)
        relation = torch.as_tensor(np.concatenate(relations), device=self.device)
        tail = torch.as_tensor(np.concatenate(tails), device=self.device)
        num_entities = num_graph_nodes + self.num_symptoms
        scorer_class = {"TransE": TransEScorer, "ComplEx": ComplExScorer}[self.model_name]
        self.scorer = scorer_class(num_entities, num_relation_types + 1, self.embedding_dim).to(self.device)
        optimizer = torch.optim.Adam(self.scorer.parameters(), lr=self.learning_rate)
        num_triples = head.shape[0]
        self.training_loss_trace = []
        for _ in range(self.num_epochs):
            order = torch.as_tensor(generator.permutation(num_triples), device=self.device)
            epoch_loss = 0.0
            for start in range(0, num_triples, self.batch_size):
                batch = order[start : start + self.batch_size]
                batch_head, batch_relation, batch_tail = head[batch], relation[batch], tail[batch]
                positive = self.scorer(batch_head, batch_relation, batch_tail)
                repeated_head = batch_head.repeat(self.negatives_per_positive)
                repeated_relation = batch_relation.repeat(self.negatives_per_positive)
                repeated_tail = batch_tail.repeat(self.negatives_per_positive)
                corrupt_tail = torch.rand(repeated_head.shape[0], device=self.device) < 0.5
                random_entities = torch.randint(0, num_entities, (repeated_head.shape[0],), device=self.device)
                negative_head = torch.where(corrupt_tail, repeated_head, random_entities)
                negative_tail = torch.where(corrupt_tail, random_entities, repeated_tail)
                negative = self.scorer(negative_head, repeated_relation, negative_tail)
                loss = torch.relu(self.margin - positive.repeat(self.negatives_per_positive) + negative).mean()
                if self.model_name == "ComplEx":
                    loss = loss + 1e-5 * sum((parameter ** 2).sum() for parameter in self.scorer.parameters()) / num_entities
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                self.scorer.normalize_entities()
                epoch_loss += loss.item() * len(batch)
            self.training_loss_trace.append(epoch_loss / num_triples)
        return self

    def predict(self, perturbation_seeds: list[np.ndarray]) -> np.ndarray:
        if self.scorer is None:
            raise RuntimeError("fit before predict")
        self.scorer.eval()
        symptom_nodes = torch.arange(self.num_graph_nodes, self.num_graph_nodes + self.num_symptoms, device=self.device)
        relation = torch.full((self.num_symptoms,), self.induces_relation_index, device=self.device, dtype=torch.long)
        scores = np.zeros((len(perturbation_seeds), self.num_symptoms))
        with torch.no_grad():
            for row, seeds in enumerate(perturbation_seeds):
                if len(seeds) == 0:
                    continue
                per_seed = []
                for seed in seeds:
                    head = torch.full((self.num_symptoms,), int(seed), device=self.device, dtype=torch.long)
                    per_seed.append(self.scorer(head, relation, symptom_nodes).cpu().numpy())
                scores[row] = np.mean(per_seed, axis=0)
        return scores
