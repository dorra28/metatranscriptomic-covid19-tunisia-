#!/usr/bin/env bash
#
# 00_run_pipeline.sh
#
# Master script — runs the full metatranscriptomic pipeline (preprocessing,
# virus analysis, host transcriptional response, microbiome profiling)
# for every sample listed in the sample sheet.
#
# Usage:
#   bash scripts/00_run_pipeline.sh config/samples.tsv [output_dir] [threads]
#
# Requires the `covid19-metatx` conda environment to be active
# (see environment.yml).

set -euo pipefail

SAMPLES_TSV="${1:?Usage: 00_run_pipeline.sh <samples.tsv> [output_dir] [threads]}"
OUTDIR="${2:-results}"
THREADS="${3:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"/{01_qc,02_virus,03_host,04_deg,05_functional,06_microbiome,07_diversity}

echo "=================================================="
echo " Metatranscriptomic COVID-19 pipeline"
echo " Samples sheet : $SAMPLES_TSV"
echo " Output dir    : $OUTDIR"
echo " Threads       : $THREADS"
echo "=================================================="

# Skip header line of the TSV
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id group r1 r2 ct rin; do

    echo ""
    echo ">>> Processing sample: $sample_id ($group)"

    # A. Preprocessing: QC + trimming
    bash "$SCRIPT_DIR/01_qc_trimming.sh" "$sample_id" "$r1" "$r2" "$OUTDIR/01_qc" "$THREADS"

    TRIMMED_R1="$OUTDIR/01_qc/${sample_id}_trimmed_R1.fastq.gz"
    TRIMMED_R2="$OUTDIR/01_qc/${sample_id}_trimmed_R2.fastq.gz"

    # B. Virus analysis
    bash "$SCRIPT_DIR/02_virus_analysis.sh" "$sample_id" "$TRIMMED_R1" "$TRIMMED_R2" "$OUTDIR/02_virus" "$THREADS"

    # C. Host mapping + quantification (per-sample; DESeq2 run separately across all samples)
    bash "$SCRIPT_DIR/03_host_mapping_quant.sh" "$sample_id" "$TRIMMED_R1" "$TRIMMED_R2" "$OUTDIR/03_host" "$THREADS"

    # D. Microbiome profiling
    bash "$SCRIPT_DIR/06_microbiome_analysis.sh" "$sample_id" "$TRIMMED_R1" "$TRIMMED_R2" "$OUTDIR/06_microbiome" "$THREADS"

done

echo ""
echo ">>> Per-sample processing complete."
echo ">>> Running cohort-level analyses (DESeq2, functional annotation, diversity)..."

# Cohort-level R analyses (run once, across all samples)
Rscript "$SCRIPT_DIR/04_deseq2_analysis.R" \
    --counts_dir "$OUTDIR/03_host" \
    --samples "$SAMPLES_TSV" \
    --outdir "$OUTDIR/04_deg"

Rscript "$SCRIPT_DIR/05_functional_annotation.R" \
    --deg_table "$OUTDIR/04_deg/DEG_results.tsv" \
    --outdir "$OUTDIR/05_functional"

Rscript "$SCRIPT_DIR/07_diversity_analysis.R" \
    --abundance_dir "$OUTDIR/06_microbiome" \
    --samples "$SAMPLES_TSV" \
    --outdir "$OUTDIR/07_diversity"

echo ""
echo "=================================================="
echo " Pipeline finished. Results in: $OUTDIR"
echo "=================================================="
