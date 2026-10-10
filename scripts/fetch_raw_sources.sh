#!/bin/bash
# Fetch every raw input that has no fetch script of its own, into data/raw/, and check each against the size and
# SHA-256 that docs/data_sources.md records. Written 10 October 2026, when a new container had to rebuild
# data/processed from nothing and seven of these inputs had a host in the documents but no URL.
#
# Each file is skipped when it already exists with a non-zero size, so the script is resumable. A recorded SHA-256 that
# does not match is reported as DRIFT and the file is kept: several sources are rolling (docs/data_sources.md), so a
# drift is a fact to record against the rebuild, not a reason to stop it. The report goes to
# data/raw/fetch_raw_sources_report.tsv with one row per file: path, bytes, sha256, recorded sha256 or "-", status.
#
# Sources with their own script are not fetched here: experiments/fetch_chembl_drug_targets.py,
# fetch_node_descriptor_sources.py (UniProt proteome, GO, hp-base.owl, ChEBI), fetch_reactome_pathways.py,
# fetch_plasma_binder_annotations.py, fetch_onsides_identifier_bridge.py, fetch_fraction_unbound_database.py,
# fetch_pubtator_relations.py, fetch_ctd_chemical_disease.py and
# python -m mechanistic_pathway_learning.evidence.gene_identifier_map (HGNC).
#
# Usage: bash scripts/fetch_raw_sources.sh [--with-cellxgene]
#   --with-cellxgene also fetches the four adult human brain atlas dissections (1.24 GB), which only the dopaminergic
#   cell class needs (experiments/build_dopaminergic_expression.py).
set -uo pipefail
repository_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repository_root"
raw_directory="data/raw"
report_path="$raw_directory/fetch_raw_sources_report.tsv"
mkdir -p "$raw_directory"
printf 'path\tbytes\tsha256\trecorded_sha256\tstatus\n' > "$report_path"
with_cellxgene="no"
if [ "${1:-}" = "--with-cellxgene" ]; then
  with_cellxgene="yes"
fi

# The Human Phenotype Ontology release docs/data_sources.md pins. The README's commands use releases/latest, which
# moves; the tag does not.
hpo_release="v2026-09-01"
hpo_release_url="https://github.com/obophenotype/human-phenotype-ontology/releases/download/$hpo_release"

fetch() {
  # fetch DESTINATION URL [RECORDED_SHA256]
  local destination_path="$1" source_url="$2" recorded_sha256="${3:--}" status="fetched"
  mkdir -p "$(dirname "$destination_path")"
  if [ -s "$destination_path" ]; then
    status="present"
  elif ! curl --fail --silent --show-error --location --retry 3 --retry-delay 5 --output "$destination_path.partial" "$source_url"; then
    rm -f "$destination_path.partial"
    printf '%s\t0\t-\t%s\tFAILED\n' "$destination_path" "$recorded_sha256" >> "$report_path"
    echo "FAILED $destination_path from $source_url" >&2
    return 1
  else
    mv "$destination_path.partial" "$destination_path"
  fi
  local measured_sha256 measured_bytes
  measured_sha256="$(sha256sum "$destination_path" | cut -d' ' -f1)"
  measured_bytes="$(stat -c %s "$destination_path")"
  if [ "$recorded_sha256" != "-" ] && [ "$recorded_sha256" != "$measured_sha256" ]; then
    status="DRIFT"
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$destination_path" "$measured_bytes" "$measured_sha256" "$recorded_sha256" "$status" >> "$report_path"
  echo "$status $destination_path ($measured_bytes bytes)"
}

# Human Phenotype Ontology: the ontology, the gene annotations, and the two files the evidence assembler reads when
# present (phenotype.hpoa for frequency and provenance, genes_to_disease.txt for the disease clusters). The assembler
# treats the last two as optional, so a missing one changes grades without an error.
for hpo_file in hp.obo genes_to_phenotype.txt phenotype.hpoa genes_to_disease.txt; do
  fetch "$raw_directory/hpo/$hpo_file" "$hpo_release_url/$hpo_file"
done

# Human-GEM 2.0.1 by its tag (commit 93e8d292a0cb). The README used to fetch the main branch, which moved to 2.1.1 on
# 7 October 2026: a build from main on 10 October gave 12,748 reactions and 8,374 metabolites against the 12,877 and
# 8,460 of data/releases/v0.4/graph_full, so the branch is not a pin.
human_gem_tag="v2.0.1"
for human_gem_file in Human-GEM.xml Human-GEM.yml genes.tsv metabolites.tsv reactions.tsv; do
  fetch "$raw_directory/Human-GEM/model/$human_gem_file" "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/$human_gem_tag/model/$human_gem_file"
done

# SIDER 4.1.
for sider_file in meddra_all_se.tsv.gz meddra_freq.tsv.gz meddra_all_indications.tsv.gz drug_names.tsv drug_atc.tsv; do
  fetch "$raw_directory/sider_4.1/$sider_file" "https://sideeffects.embl.de/media/download/$sider_file"
done

# OmniPath web service: signalling, CollecTRI regulons and small-molecule interactions. A live service with no
# release tag, so the row counts of docs/data_sources.md (126,124, 64,515 and 6,553) are the only check.
fetch "$raw_directory/omnipath/omnipath_interactions.tsv" "https://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,ligrecextra&genesymbols=1&fields=sources,references,type,curation_effort,n_references&organisms=9606&format=tsv"
fetch "$raw_directory/omnipath/collectri_interactions.tsv" "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=sources,references,n_references&organisms=9606&format=tsv"
fetch "$raw_directory/omnipath/small_molecule_protein.tsv" "https://omnipathdb.org/interactions?types=small_molecule_protein&genesymbols=1&fields=sources,references,type&organisms=9606&format=tsv"

# GTEx v10 median gene TPM (8,846,936 bytes recorded in docs/alternative_sources.md).
fetch "$raw_directory/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"

# Human Protein Atlas: brain regions and single-nucleus cluster types.
fetch "$raw_directory/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip" "https://www.proteinatlas.org/download/tsv/rna_brain_region_hpa.tsv.zip" 108fb699eccc3f8611b475caa7968ad67f10ee47c4cc88db16d2fcaa0214a0f3
fetch "$raw_directory/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip" "https://www.proteinatlas.org/download/tsv/rna_single_nuclei_cluster_type.tsv.zip" a287c1abab593c9d6802349a82f118efc0d4aea411be7d897b8feb0b31b02f2a

# Orphadata product 6: gene-disease association types.
fetch "$raw_directory/orphadata/en_product6.xml" "https://www.orphadata.com/data/xml/en_product6.xml"

# OnSIDES v3.1.1 release asset (84,862,297 bytes recorded in docs/alternative_sources.md). The builds read parquet
# tables under data/raw/onsides/v3.1.1/parquet; scripts/convert_onsides_to_parquet.py makes them from this archive.
fetch "$raw_directory/onsides/v3.1.1/onsides-v3.1.1.zip" "https://github.com/tatonetti-lab/onsides/releases/download/v3.1.1/onsides-v3.1.1.zip"

# ESM-2 esm2_t12_35M_UR50D weights and contact regression.
fetch "$raw_directory/esm/esm2_t12_35M_UR50D.pt" "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t12_35M_UR50D.pt" 7f21e80e61d16a71735163ef555d3009afb0c98da74c48e29df08606973cc55e
fetch "$raw_directory/esm/esm2_t12_35M_UR50D-contact-regression.pt" "https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t12_35M_UR50D-contact-regression.pt" 16641e05d830d0ce863dd152dbb8c2f3ddfa3c3ec2a66080152c8abad01d8585

# HGNC complete set, from the monthly archive copy that carries the pinned hash. The rolling URL the gene identifier
# map downloads from served a later file on 10 October 2026 (16,973,159 bytes against the pinned 16,963,116), and the
# map keeps a cached file that has a pin beside it, so fetching the archive copy first is what keeps the pin.
fetch "$raw_directory/hgnc/hgnc_complete_set.txt" "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-02.txt" 2b4224ea847df2fc6982f5b2a52804c5d92fbb8810b134afb636f6029452dc03
if [ ! -s "$raw_directory/hgnc/hgnc_complete_set.pin.json" ]; then
  python3 - "$raw_directory/hgnc" <<'PIN'
import hashlib, json, sys
from pathlib import Path
hgnc_path = Path(sys.argv[1]) / "hgnc_complete_set.txt"
pin = {"url": "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-02.txt",
       "last_modified": "Fri, 02 Oct 2026 13:44:59 GMT", "content_length": hgnc_path.stat().st_size,
       "sha256": hashlib.sha256(hgnc_path.read_bytes()).hexdigest(), "row_count": sum(1 for _ in hgnc_path.open()) - 1}
(Path(sys.argv[1]) / "hgnc_complete_set.pin.json").write_text(json.dumps(pin, indent=2))
PIN
fi

# UniProt reviewed human proteins whose cofactor annotation names Mn(2+) (ChEBI 29035), for the manganese variant
# (experiments/build_manganese_graph_variant.py). The old container held this file with no script behind it; the
# query below returned 281 entries on release 2026_03, which is the 93 requiring manganese alone plus the 188
# accepting it among other metals that the variant's summary records.
fetch "$raw_directory/uniprot/uniprot_human_manganese_cofactor.tsv.gz" "https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)%20AND%20(reviewed:true)%20AND%20(cc_cofactor_chebi:%22CHEBI:29035%22)&fields=accession,gene_primary,cc_cofactor,cc_subcellular_location&format=tsv&compressed=true"

if [ "$with_cellxgene" = "yes" ]; then
  # Adult human brain atlas dissections holding cluster 395 (Siletti et al. 2023; CELLxGENE collection
  # 283d65eb-dd53-496d-adb7-7570c7caa443). The identifiers docs/data_sources.md lists are dataset version ids, which
  # name one immutable file each, and the file is served under that id. docs/data_sources.md records the first 16 hex
  # digits of each SHA-256, so those are compared here and printed, since fetch() compares whole digests only.
  declare -A dissection_versions=(
    [dissection_sn_rn]="0beb5c17-8951-4e94-a5d9-9b8e13784395:96b6f96d706317a6"
    [dissection_sn]="2de1ee2f-d82a-4bea-9288-bae4ea063138:b879112e04a3eedb"
    [dissection_pag_dr]="5b5716d6-8900-478c-824d-3cc7bf2c99ce:4460667030694cd4"
    [dissection_pag]="0405d1d6-dd88-49d6-9aac-fdb6699612cd:7be3ae310f3e70df"
  )
  for dissection_name in dissection_sn_rn dissection_sn dissection_pag_dr dissection_pag; do
    version_identifier="${dissection_versions[$dissection_name]%%:*}"
    recorded_prefix="${dissection_versions[$dissection_name]##*:}"
    dissection_path="$raw_directory/siletti_cellxgene/$dissection_name.h5ad"
    fetch "$dissection_path" "https://datasets.cellxgene.cziscience.com/$version_identifier.h5ad" || continue
    measured_prefix="$(sha256sum "$dissection_path" | cut -c1-16)"
    if [ "$measured_prefix" = "$recorded_prefix" ]; then
      echo "  sha256 prefix $measured_prefix matches the recorded one"
    else
      echo "  DRIFT: sha256 prefix $measured_prefix, recorded $recorded_prefix" >&2
      printf '%s\t-\t%s\t%s\tDRIFT\n' "$dissection_path" "$measured_prefix" "$recorded_prefix" >> "$report_path"
    fi
  done
fi

echo "report: $report_path"
if grep -q -P '\tFAILED$' "$report_path"; then
  exit 1
fi
