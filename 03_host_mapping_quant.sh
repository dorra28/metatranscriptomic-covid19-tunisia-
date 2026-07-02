#!/usr/bin/env bash
#
# 03_host_mapping_quant.sh — Host transcriptional response, mapping + quantification
# (Figure 1, panel C, steps 1-3)
#
# Maps trimmed reads to the human reference genome (GRCh38) and quantifies
# gene-level counts with FeatureCounts. Per-sample count files are combined
# and passed to DESeq2 in 04_deseq2_analysis.R.
#
# Tools / versions: BWA-MEM 0.7.17, Samtools, FeatureCounts (Subread) 2.0.1
#
# Usage:
#   bash 03_host_mapping_quant.sh <sample_id> <trimmed_R1> <trimmed_R2> <outdir> [threads]
#
# Required reference files:
#   references/GRCh38.fasta      — human reference genome, BWA-indexed
#   references/GRCh38.gtf        — gene annotation (GENCODE/Ensembl), for FeatureCounts

set -euo pipefail

SAMPLE_ID="${1:?sample_id required}"
R1="${2:?trimmed R1 required}"
R2="${3:?trimmed R2 required}"
OUTDIR="${4:?outdir required}"
THREADS="${5:-8}"

REF="${REF_GRCH38:-references/GRCh38.fasta}"
GTF="${GTF_GRCH38:-references/GRCh38.gtf}"

mkdir -p "$OUTDIR"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "[03_host_mapping_quant] Indexing GRCh38 reference (this can take a while)..."
    bwa index "$REF"
fi

echo "[03_host_mapping_quant] $SAMPLE_ID — BWA-MEM mapping to GRCh38"
bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE_ID}.grch38.sorted.bam" -
samtools index "$OUTDIR/${SAMPLE_ID}.grch38.sorted.bam"

echo "[03_host_mapping_quant] $SAMPLE_ID — FeatureCounts gene-level quantification"
featureCounts \
    -T "$THREADS" \
    -p --countReadPairs \
    -a "$GTF" \
    -o "$OUTDIR/${SAMPLE_ID}.featureCounts.txt" \
    "$OUTDIR/${SAMPLE_ID}.grch38.sorted.bam"

# Keep a simplified two-column (gene_id, count) file for easy merging in R
tail -n +3 "$OUTDIR/${SAMPLE_ID}.featureCounts.txt" \
    | awk -v OFS='\t' '{print $1, $NF}' \
    > "$OUTDIR/${SAMPLE_ID}.counts.tsv"

echo "[03_host_mapping_quant] $SAMPLE_ID — done"
