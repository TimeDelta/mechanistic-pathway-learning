"""ESM-2 (35M parameters, 480 dimensions) embeddings of the reviewed human proteome, for the protein descriptors
(mechanistic_pathway_learning/graph/protein_descriptors.py; design section 5.2).

Each protein is embedded as the mean over its residues of the last layer (layer 12); a sequence longer than the model's
1,022-residue context is cut into consecutive windows of at most 1,022 residues and the window means are averaged with
weights equal to their lengths. Proteins of graph genes come first, so a partial run already covers the slice.

The job is idempotent and resumable: proteins are processed in batches of --batch-proteins, each written to
<output-dir>/batch<j>.npz (accessions and embeddings) and skipped when present; --shard-index and --num-shards
split the work between processes. The weights are the pinned files of data/raw/esm (loaded with weights_only and only
argparse.Namespace allowlisted, the one non-tensor object the checkpoint holds).

Usage:
  OMP_NUM_THREADS=1 python experiments/compute_protein_embeddings.py --shard-index 0 --num-shards 2
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

MODEL_NAME = "esm2_t12_35M_UR50D"
REPRESENTATION_LAYER = 12
MAXIMUM_WINDOW_RESIDUES = 1022


def load_esm_model(weights_directory: Path):
    import argparse as argparse_module

    import esm

    with torch.serialization.safe_globals([argparse_module.Namespace]):
        model_data = torch.load(weights_directory / f"{MODEL_NAME}.pt", map_location="cpu", weights_only=True)
        regression_data = torch.load(weights_directory / f"{MODEL_NAME}-contact-regression.pt", map_location="cpu", weights_only=True)
    model, alphabet = esm.pretrained.load_model_and_alphabet_core(MODEL_NAME, model_data, regression_data)
    model.eval()
    return model, alphabet


def embed_sequence(model, batch_converter, accession: str, sequence: str) -> np.ndarray:
    window_means, window_lengths = [], []
    with torch.no_grad():
        for start in range(0, len(sequence), MAXIMUM_WINDOW_RESIDUES):
            window = sequence[start:start + MAXIMUM_WINDOW_RESIDUES]
            _, _, tokens = batch_converter([(accession, window)])
            residues = model(tokens, repr_layers=[REPRESENTATION_LAYER])["representations"][REPRESENTATION_LAYER][0, 1:len(window) + 1]
            window_means.append(residues.mean(dim=0).numpy())
            window_lengths.append(len(window))
    return np.average(np.stack(window_means), axis=0, weights=np.asarray(window_lengths, dtype=float)).astype(np.float32)


def ordered_proteins(uniprot_table: pd.DataFrame, graph_gene_symbols: set[str]) -> pd.DataFrame:
    """Graph genes first, then the rest, each by accession, so the batch layout is fixed across runs."""
    table = uniprot_table.assign(in_graph=uniprot_table["Gene Names (primary)"].isin(graph_gene_symbols))
    return table.sort_values(["in_graph", "Entry"], ascending=[False, True]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--graph-nodes", type=Path, default=Path("data/processed/graph_full/nodes.parquet"))
    parser.add_argument("--weights-directory", type=Path, default=Path("data/raw/esm"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/protein_embeddings"))
    parser.add_argument("--batch-proteins", type=int, default=200)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--num-shards", type=int, default=1)
    arguments = parser.parse_args()
    torch.set_num_threads(1)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    nodes = pd.read_parquet(arguments.graph_nodes)
    proteins = ordered_proteins(pd.read_csv(arguments.uniprot_table, sep="\t"), set(nodes.gene_symbol.dropna()))
    batches = [proteins.iloc[start:start + arguments.batch_proteins] for start in range(0, len(proteins), arguments.batch_proteins)]
    own_batches = [(index, batch) for index, batch in enumerate(batches) if index % arguments.num_shards == arguments.shard_index]
    model, alphabet = None, None
    started = time.time()
    for done, (batch_index, batch) in enumerate(own_batches):
        output_path = arguments.output_dir / f"batch{batch_index:04d}.npz"
        if output_path.exists():
            continue
        if model is None:
            model, alphabet = load_esm_model(arguments.weights_directory)
            batch_converter = alphabet.get_batch_converter()
        embeddings = np.stack([embed_sequence(model, batch_converter, accession, sequence) for accession, sequence in zip(batch.Entry, batch.Sequence)])
        partial_path = arguments.output_dir / f"batch{batch_index:04d}.part.npz"
        np.savez(partial_path, accessions=batch.Entry.to_numpy(dtype=str), embeddings=embeddings)
        partial_path.replace(output_path)
        print(f"batch {batch_index} ({done + 1} of {len(own_batches)} in shard {arguments.shard_index}) done, {time.time() - started:.0f}s", flush=True)
    finished = sum((arguments.output_dir / f"batch{index:04d}.npz").exists() for index in range(len(batches)))
    print(f"{finished} of {len(batches)} batches present", flush=True)


if __name__ == "__main__":
    main()
