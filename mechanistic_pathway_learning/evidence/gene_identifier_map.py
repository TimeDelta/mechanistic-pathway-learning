"""Gene symbol <-> NCBI Gene id <-> Ensembl gene id for every gene node of the physiology graph (design section 4.2, E3).

Why: PubTator3 names genes by NCBI Gene id and the graph names them by symbol (metabolic genes from Human-GEM also
carry an Ensembl id; signaling and transcription genes from OmniPath and CollecTRI carry only the symbol). The
literature retrieval needs the NCBI id of every graph gene, so the join runs on a stable integer rather than on a
symbol that may have been renamed since the paper was annotated.

Precedence, one source per node, recorded in mapping_source:
  1. human_gem_ensembl       Human-GEM genes.tsv row whose Ensembl id equals the node's Ensembl id
  2. human_gem_symbol        Human-GEM genes.tsv row whose symbol equals the node's symbol
  3. hgnc_symbol             HGNC complete set, current approved symbol
  4. hgnc_previous_symbol    HGNC previous symbol, only when exactly one current gene carried it
  5. hgnc_alias_symbol       HGNC alias symbol, only when exactly one current gene carries it
  6. unmapped                no NCBI id; the node gets no literature reports
Human-GEM comes first because the graph's metabolic layer was built from it, so its ids are the ones the
reconstruction itself asserts; HGNC is the authority for the rest. A previous or alias symbol shared by several
genes is ambiguous and is left unmapped rather than guessed.

The HGNC complete set is downloaded once to data/raw/hgnc (gitignored) with a pin file recording the URL, size,
Last-Modified header, SHA-256 and fetch date.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

HGNC_COMPLETE_SET_URLS: tuple[str, ...] = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
)
HGNC_FILE_NAME = "hgnc_complete_set.txt"
HGNC_PIN_FILE_NAME = "hgnc_complete_set.pin.json"
MAPPING_SOURCES: tuple[str, ...] = ("human_gem_ensembl", "human_gem_symbol", "hgnc_symbol", "hgnc_previous_symbol", "hgnc_alias_symbol", "unmapped")
GENE_MAP_COLUMNS: tuple[str, ...] = ("gene_symbol", "ensembl_gene_id", "ncbi_gene_id", "mapping_source", "mapped_symbol")


@dataclass(frozen=True)
class GeneIdentifierRow:
    gene_symbol: str  # the graph node's symbol
    ensembl_gene_id: str | None
    ncbi_gene_id: int | None
    mapping_source: str
    mapped_symbol: str | None  # the symbol of the source row that supplied the ids (differs from gene_symbol for previous or alias matches)


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_hgnc_complete_set(cache_directory: Path, urls: tuple[str, ...] = HGNC_COMPLETE_SET_URLS) -> Path:
    """The cached HGNC file, downloaded from the first URL that answers when absent; writes the pin file."""
    cache_directory = Path(cache_directory)
    cache_directory.mkdir(parents=True, exist_ok=True)
    target_path = cache_directory / HGNC_FILE_NAME
    pin_path = cache_directory / HGNC_PIN_FILE_NAME
    if target_path.exists() and pin_path.exists():
        return target_path
    last_error: Exception | None = None
    for url in urls:
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "mechanistic-pathway-learning/0.1"})
            with urllib.request.urlopen(request, timeout=120) as response:
                headers = {key.lower(): value for key, value in response.headers.items()}
                temporary_path = target_path.with_suffix(".txt.part")
                with open(temporary_path, "wb") as handle:
                    for chunk in iter(lambda: response.read(1 << 20), b""):
                        handle.write(chunk)
                temporary_path.replace(target_path)
            pin = {
                "url": url,
                "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "last_modified": headers.get("last-modified"),
                "content_length": target_path.stat().st_size,
                "sha256": sha256_of_file(target_path),
                "row_count": sum(1 for _ in open(target_path, encoding="utf-8")) - 1,
            }
            pin_path.write_text(json.dumps(pin, indent=2))
            return target_path
        except Exception as error:  # noqa: BLE001 - try the mirror
            last_error = error
    raise RuntimeError(f"could not download the HGNC complete set from {urls}") from last_error


def write_hgnc_pin_for_existing_file(cache_directory: Path, url: str, last_modified: str | None) -> Path:
    """Pin a file that was downloaded outside this module (the headers are passed in); idempotent."""
    cache_directory = Path(cache_directory)
    target_path = cache_directory / HGNC_FILE_NAME
    pin_path = cache_directory / HGNC_PIN_FILE_NAME
    if not pin_path.exists():
        pin = {
            "url": url,
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(target_path.stat().st_mtime)),
            "last_modified": last_modified,
            "content_length": target_path.stat().st_size,
            "sha256": sha256_of_file(target_path),
            "row_count": sum(1 for _ in open(target_path, encoding="utf-8")) - 1,
        }
        pin_path.write_text(json.dumps(pin, indent=2))
    return pin_path


def load_human_gem_gene_table(genes_tsv_path: Path) -> pd.DataFrame:
    """Columns ensembl_gene_id, gene_symbol, ncbi_gene_id (Int64) from Human-GEM model/genes.tsv."""
    table = pd.read_csv(genes_tsv_path, sep="\t", dtype=str, quoting=csv.QUOTE_ALL, keep_default_na=False)
    return pd.DataFrame(
        {
            "ensembl_gene_id": table["genes"].str.strip(),
            "gene_symbol": table["geneSymbols"].str.strip(),
            "ncbi_gene_id": pd.to_numeric(table["geneEntrezID"].str.strip().replace("", None), errors="coerce").astype("Int64"),
        }
    )


def load_hgnc_table(hgnc_path: Path) -> pd.DataFrame:
    """Approved HGNC genes with columns symbol, ncbi_gene_id (Int64), ensembl_gene_id, prev_symbols, alias_symbols (lists)."""
    table = pd.read_csv(hgnc_path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    approved = table[table["status"] == "Approved"] if "status" in table.columns else table
    return pd.DataFrame(
        {
            "symbol": approved["symbol"].str.strip().values,
            "ncbi_gene_id": pd.to_numeric(approved["entrez_id"].str.strip().replace("", None), errors="coerce").astype("Int64").values,
            "ensembl_gene_id": approved["ensembl_gene_id"].str.strip().replace("", None).values,
            "prev_symbols": [[entry.strip() for entry in cell.split("|") if entry.strip()] for cell in approved["prev_symbol"]],
            "alias_symbols": [[entry.strip() for entry in cell.split("|") if entry.strip()] for cell in approved["alias_symbol"]],
        }
    )


def _unique_symbol_index(hgnc: pd.DataFrame, list_column: str) -> dict[str, int]:
    """Symbol -> row position, only for symbols that exactly one HGNC gene lists in that column."""
    positions: dict[str, list[int]] = defaultdict(list)
    for position, symbols in enumerate(hgnc[list_column]):
        for symbol in symbols:
            positions[symbol].append(position)
    return {symbol: rows[0] for symbol, rows in positions.items() if len(rows) == 1}


def map_gene_nodes(gene_nodes: pd.DataFrame, human_gem: pd.DataFrame, hgnc: pd.DataFrame) -> pd.DataFrame:
    """One row per gene node with the ids resolved by the precedence in the module docstring.

    gene_nodes needs columns gene_symbol and ensembl_gene_id (null allowed); human_gem and hgnc are the tables of the
    two loaders above. Matching is on exact symbols: the graph's symbols are HGNC symbols by construction, so
    case-folding would only invite false matches between distinct genes.
    """
    human_gem_by_ensembl = {row.ensembl_gene_id: row for row in human_gem.itertuples(index=False) if row.ensembl_gene_id}
    human_gem_by_symbol: dict[str, object] = {}
    for row in human_gem.itertuples(index=False):
        human_gem_by_symbol.setdefault(row.gene_symbol, row)
    hgnc_by_symbol = {row.symbol: row for row in hgnc.itertuples(index=False)}
    hgnc_rows = list(hgnc.itertuples(index=False))
    previous_index = _unique_symbol_index(hgnc, "prev_symbols")
    alias_index = _unique_symbol_index(hgnc, "alias_symbols")

    def integer_or_none(value) -> int | None:
        return None if value is None or pd.isna(value) else int(value)

    def text_or_none(value) -> str | None:
        return None if value is None or (isinstance(value, float) and pd.isna(value)) or value == "" else str(value)

    mapped: list[GeneIdentifierRow] = []
    for node in gene_nodes.itertuples(index=False):
        symbol = str(node.gene_symbol)
        node_ensembl = text_or_none(node.ensembl_gene_id)
        source_row = None
        mapping_source = "unmapped"
        if node_ensembl is not None and node_ensembl in human_gem_by_ensembl:
            source_row, mapping_source = human_gem_by_ensembl[node_ensembl], "human_gem_ensembl"
        elif symbol in human_gem_by_symbol:
            source_row, mapping_source = human_gem_by_symbol[symbol], "human_gem_symbol"
        elif symbol in hgnc_by_symbol:
            source_row, mapping_source = hgnc_by_symbol[symbol], "hgnc_symbol"
        elif symbol in previous_index:
            source_row, mapping_source = hgnc_rows[previous_index[symbol]], "hgnc_previous_symbol"
        elif symbol in alias_index:
            source_row, mapping_source = hgnc_rows[alias_index[symbol]], "hgnc_alias_symbol"
        if source_row is None:
            mapped.append(GeneIdentifierRow(symbol, node_ensembl, None, "unmapped", None))
            continue
        ncbi_gene_id = integer_or_none(source_row.ncbi_gene_id)
        if ncbi_gene_id is None:
            # the matched row has no NCBI id (some Human-GEM rows); fall back to HGNC by the same symbol when it helps
            hgnc_row = hgnc_by_symbol.get(getattr(source_row, "gene_symbol", getattr(source_row, "symbol", symbol)))
            if hgnc_row is not None and integer_or_none(hgnc_row.ncbi_gene_id) is not None:
                source_row, mapping_source, ncbi_gene_id = hgnc_row, "hgnc_symbol", integer_or_none(hgnc_row.ncbi_gene_id)
        source_symbol = getattr(source_row, "gene_symbol", None) or getattr(source_row, "symbol", None)
        ensembl_gene_id = node_ensembl or text_or_none(source_row.ensembl_gene_id)
        mapped.append(GeneIdentifierRow(symbol, ensembl_gene_id, ncbi_gene_id, mapping_source if ncbi_gene_id is not None else "unmapped", source_symbol))
    table = pd.DataFrame([row.__dict__ for row in mapped], columns=list(GENE_MAP_COLUMNS))
    table["ncbi_gene_id"] = table["ncbi_gene_id"].astype("Int64")
    return table


def ncbi_gene_id_to_symbol(gene_map: pd.DataFrame) -> dict[int, str]:
    """NCBI id -> graph symbol; when two graph nodes resolve to one NCBI id (a renamed gene present under both
    names) the first in table order wins and the duplicate is reported by the caller."""
    lookup: dict[int, str] = {}
    for row in gene_map.itertuples(index=False):
        if pd.notna(row.ncbi_gene_id) and int(row.ncbi_gene_id) not in lookup:
            lookup[int(row.ncbi_gene_id)] = row.gene_symbol
    return lookup


def build_gene_identifier_map(graph_directory: Path, human_gem_genes_path: Path, hgnc_directory: Path, output_path: Path | None) -> pd.DataFrame:
    nodes = pd.read_parquet(Path(graph_directory) / "nodes.parquet", columns=["node_id", "node_type", "gene_symbol", "ensembl_gene_id"])
    gene_nodes = nodes[nodes.node_type == "gene"][["gene_symbol", "ensembl_gene_id"]].reset_index(drop=True)
    human_gem = load_human_gem_gene_table(human_gem_genes_path)
    hgnc = load_hgnc_table(download_hgnc_complete_set(hgnc_directory))
    gene_map = map_gene_nodes(gene_nodes, human_gem, hgnc)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        gene_map.to_parquet(output_path, index=False)
    return gene_map


def gene_map_summary(gene_map: pd.DataFrame) -> dict:
    ncbi_known = gene_map[gene_map.ncbi_gene_id.notna()]
    duplicated_ncbi = ncbi_known[ncbi_known.ncbi_gene_id.duplicated(keep=False)]
    return {
        "gene_nodes": int(len(gene_map)),
        "with_ncbi_gene_id": int(len(ncbi_known)),
        "without_ncbi_gene_id": int(len(gene_map) - len(ncbi_known)),
        "with_ensembl_gene_id": int(gene_map.ensembl_gene_id.notna().sum()),
        "by_mapping_source": {source: int(count) for source, count in gene_map.mapping_source.value_counts().sort_index().items()},
        "ncbi_ids_shared_by_several_nodes": int(duplicated_ncbi.ncbi_gene_id.nunique()),
        "nodes_on_shared_ncbi_ids": sorted(duplicated_ncbi.gene_symbol.tolist()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full"))
    parser.add_argument("--human-gem-genes", type=Path, default=Path("data/raw/Human-GEM/model/genes.tsv"))
    parser.add_argument("--hgnc-dir", type=Path, default=Path("data/raw/hgnc"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/literature/gene_identifier_map.parquet"))
    arguments = parser.parse_args()
    gene_map = build_gene_identifier_map(arguments.graph_dir, arguments.human_gem_genes, arguments.hgnc_dir, arguments.output)
    summary = gene_map_summary(gene_map)
    summary_path = arguments.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
