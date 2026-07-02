#!/usr/bin/env bash
#
# 06_microbiome_analysis.sh — Microbiome profiling (Figure 1, panel D)
#
# Removes host/viral reads, filters chimeras, and performs taxonomic
# profiling of the remaining (microbial) reads.
#
# Tools / versions: Bowtie2 2.4.5, VSEARCH 2.22.0, MetaPhlAn 4.0.6
#
# Usage:
#   bash 06_microbiome_analysis.sh <sample_id> <trimmed_R1> <trimmed_R2> <outdir> [threads]
#
# Required reference files:
#   references/host_viral_bowtie2_index/  — Bowtie2 index built from GRCh38 + NC_045512.2
#                                            (combined host + viral reference, for read removal)
#   MetaPhlAn database is fetched automatically on first run via
#   `metaphlan --install`, or point to an existing DB with --bowtie2db.

set -euo pipefail

SAMPLE_ID="${1:?sample_id required}"
R1="${2:?trimmed R1 required}"
R2="${3:?trimmed R2 required}"
OUTDIR="${4:?outdir required}"
THREADS="${5:-8}"

HOST_VIRAL_INDEX="${HOST_VIRAL_BOWTIE2_INDEX:-references/host_viral_bowtie2_index/combined}"
METAPHLAN_DB="${METAPHLAN_DB:-}"   # optional: path to pre-installed MetaPhlAn db

mkdir -p "$OUTDIR/$SAMPLE_ID"
cd "$OUTDIR/$SAMPLE_ID"

echo "[06_microbiome_analysis] $SAMPLE_ID — Bowtie2 removal of host/viral reads"
bowtie2 \
    -x "$HOST_VIRAL_INDEX" \
    -1 "$R1" -2 "$R2" \
    -p "$THREADS" \
    --un-conc-gz "${SAMPLE_ID}_unmapped_R%.fastq.gz" \
    -S /dev/null \
    2> "${SAMPLE_ID}_bowtie2.log"

mv "${SAMPLE_ID}_unmapped_R1.fastq.gz" "${SAMPLE_ID}_microbial_R1.fastq.gz"
mv "${SAMPLE_ID}_unmapped_R2.fastq.gz" "${SAMPLE_ID}_microbial_R2.fastq.gz"

echo "[06_microbiome_analysis] $SAMPLE_ID — VSEARCH chimera filtering"
# VSEARCH chimera detection works on merged/single-end reads; merge pairs first
vsearch --fastq_mergepairs "${SAMPLE_ID}_microbial_R1.fastq.gz" \
    --reverse "${SAMPLE_ID}_microbial_R2.fastq.gz" \
    --fastqout "${SAMPLE_ID}_merged.fastq" \
    --threads "$THREADS" \
    2> "${SAMPLE_ID}_vsearch_merge.log"

vsearch --fastx_filter "${SAMPLE_ID}_merged.fastq" \
    --fastaout "${SAMPLE_ID}_merged.fasta" \
    2> "${SAMPLE_ID}_vsearch_convert.log"

vsearch --uchime_denovo "${SAMPLE_ID}_merged.fasta" \
    --nonchimeras "${SAMPLE_ID}_nonchimeric.fasta" \
    --chimeras "${SAMPLE_ID}_chimeras.fasta" \
    --threads "$THREADS" \
    2> "${SAMPLE_ID}_vsearch_chimera.log"

echo "[06_microbiome_analysis] $SAMPLE_ID — MetaPhlAn4 taxonomic profiling"
METAPHLAN_ARGS=(--input_type fasta --nproc "$THREADS")
if [[ -n "$METAPHLAN_DB" ]]; then
    METAPHLAN_ARGS+=(--bowtie2db "$METAPHLAN_DB")
fi

metaphlan "${SAMPLE_ID}_nonchimeric.fasta" \
    "${METAPHLAN_ARGS[@]}" \
    --bowtie2out "${SAMPLE_ID}.bowtie2.bz2" \
    -o "${SAMPLE_ID}_metaphlan_profile.tsv"

echo "[06_microbiome_analysis] $SAMPLE_ID — done"
