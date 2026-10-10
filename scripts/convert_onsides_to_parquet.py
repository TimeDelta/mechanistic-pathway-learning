"""Convert the OnSIDES v3.1.1 release archive into the parquet tables the label loader reads
(mechanistic_pathway_learning/evidence/load_onsides_label_events.py, load_onsides_tables).

The archive ships CSV; the loader and experiments/fetch_onsides_identifier_bridge.py read
data/raw/onsides/v3.1.1/parquet/<table>.parquet. The conversion existed only on the disk of the container that first
ran it, so it is written down here (10 October 2026). Two choices in it are not free:

- label_section keeps the literal string "NA", which is OnSIDES' value for the undivided non-US label. A reader that
  treats "NA" as missing (pandas.read_csv by default) drops every EU, UK and JP statement at the loader's group-by.
- RxNorm identifiers are text, because the vocabulary mixes RxNorm numbers with OMOP extension identifiers such as
  "OMOP997977"; label and MedDRA identifiers are integers. These are the types of the fixture in
  tests/test_onsides_label_events.py.

Usage:
  python scripts/convert_onsides_to_parquet.py        # reads data/raw/onsides/v3.1.1/onsides-v3.1.1.zip
"""
from __future__ import annotations

import argparse
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pyarrow_csv
import pyarrow.parquet as pyarrow_parquet

COLUMN_TYPES_BY_TABLE: dict[str, dict[str, pa.DataType]] = {
    "product_adverse_effect": {"product_label_id": pa.int64(), "effect_id": pa.int64(), "label_section": pa.string(), "effect_meddra_id": pa.int64(),
                               "match_method": pa.string(), "pred0": pa.float64(), "pred1": pa.float64()},
    "product_label": {"label_id": pa.int64(), "source": pa.string(), "source_product_name": pa.string(), "source_product_id": pa.string(),
                      "source_label_url": pa.string()},
    "product_to_rxnorm": {"label_id": pa.int64(), "rxnorm_product_id": pa.string()},
    "vocab_rxnorm_ingredient_to_product": {"product_id": pa.string(), "ingredient_id": pa.string()},
    "vocab_rxnorm_ingredient": {"rxnorm_id": pa.string(), "rxnorm_name": pa.string(), "rxnorm_term_type": pa.string()},
    "vocab_rxnorm_product": {"rxnorm_id": pa.string(), "rxnorm_name": pa.string(), "rxnorm_term_type": pa.string()},
    "vocab_meddra_adverse_effect": {"meddra_id": pa.int64(), "meddra_name": pa.string(), "meddra_term_type": pa.string()},
    "high_confidence": {"ingredient_id": pa.string(), "effect_meddra_id": pa.int64()},
}
# only numeric columns may read an empty field as missing; a text column keeps every value as written, "NA" included
NUMERIC_NULL_SPELLINGS: list[str] = [""]


def convert_archive(archive_path: Path, output_directory: Path) -> dict[str, int]:
    """Row count per table written to output_directory."""
    output_directory.mkdir(parents=True, exist_ok=True)
    row_counts: dict[str, int] = {}
    with zipfile.ZipFile(archive_path) as archive:
        for table_name, column_types in COLUMN_TYPES_BY_TABLE.items():
            with archive.open(f"csv/{table_name}.csv") as table_stream:
                table = pyarrow_csv.read_csv(
                    table_stream,
                    parse_options=pyarrow_csv.ParseOptions(newlines_in_values=True),
                    convert_options=pyarrow_csv.ConvertOptions(column_types=column_types, null_values=NUMERIC_NULL_SPELLINGS, strings_can_be_null=False),
                )
            pyarrow_parquet.write_table(table, output_directory / f"{table_name}.parquet")
            row_counts[table_name] = table.num_rows
    return row_counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--archive", type=Path, default=Path("data/raw/onsides/v3.1.1/onsides-v3.1.1.zip"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw/onsides/v3.1.1/parquet"))
    arguments = parser.parse_args()
    for table_name, row_count in convert_archive(arguments.archive, arguments.output_dir).items():
        print(f"{table_name}: {row_count:,} rows")


if __name__ == "__main__":
    main()
