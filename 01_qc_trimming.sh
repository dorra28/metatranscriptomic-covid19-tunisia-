#!/usr/bin/env bash
#
# 01_qc_trimming.sh — Preprocessing step (Figure 1, panel A)
#
# Quality control (FastQC + MultiQC) and adapter trimming (Cutadapt)
# of raw paired-end FASTQ files.
#
# Tools / versions: FastQC 0.12.0, MultiQC 1.14, Cutadapt 1.18
#
# Usage:
#   bash 01_qc_trimming.sh <sample_id> <R1.fastq.gz> <R2.fastq.gz> <outdir> [threads]

set -euo pipefail

SAMPLE_ID="${1:?sample_id required}"
R1="${2:?R1 fastq required}"
R2="${3:?R2 fastq required}"
OUTDIR="${4:?outdir required}"
THREADS="${5:-8}"

mkdir -p "$OUTDIR/fastqc_raw" "$OUTDIR/fastqc_trimmed"

echo "[01_qc_trimming] $SAMPLE_ID — raw FastQC"
fastqc -t "$THREADS" -o "$OUTDIR/fastqc_raw" "$R1" "$R2"

echo "[01_qc_trimming] $SAMPLE_ID — Cutadapt trimming"
# Standard Illumina TruSeq adapters — replace with your library prep kit's
# adapter sequences if different.
cutadapt \
    -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -q 20,20 \
    -m 36 \
    -j "$THREADS" \
    -o "$OUTDIR/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
    -p "$OUTDIR/${SAMPLE_ID}_trimmed_R2.fastq.gz" \
    "$R1" "$R2" \
    > "$OUTDIR/${SAMPLE_ID}_cutadapt.log" 2>&1

echo "[01_qc_trimming] $SAMPLE_ID — trimmed FastQC"
fastqc -t "$THREADS" -o "$OUTDIR/fastqc_trimmed" \
    "$OUTDIR/${SAMPLE_ID}_trimmed_R1.fastq.gz" \
    "$OUTDIR/${SAMPLE_ID}_trimmed_R2.fastq.gz"

echo "[01_qc_trimming] $SAMPLE_ID — MultiQC aggregate report"
multiqc -f -o "$OUTDIR" "$OUTDIR/fastqc_raw" "$OUTDIR/fastqc_trimmed" \
    -n "${SAMPLE_ID}_multiqc_report"

echo "[01_qc_trimming] $SAMPLE_ID — done"
