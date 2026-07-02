#!/usr/bin/env bash
#
# 02_virus_analysis.sh — SARS-CoV-2 analysis (Figure 1, panel B)
#
# Mapping to the SARS-CoV-2 reference (Wuhan-Hu-1, NC_045512.2), mapping QC,
# genome assembly, alignment, lineage/clade assignment, and variant calling.
#
# Tools / versions: BWA-MEM 0.7.17, Samtools, Qualimap 2.2.1, SPAdes 3.15.3,
#                    MAFFT 7.471, Nextclade 1.5.3, LoFreq 2.1.5, SnpEff
#
# Thresholds: minimum genome coverage 85%, minimum depth >51X for variant calling.
#
# Usage:
#   bash 02_virus_analysis.sh <sample_id> <trimmed_R1> <trimmed_R2> <outdir> [threads]
#
# Required reference files (download separately, see docs/methods_parameters.md):
#   references/NC_045512.2.fasta       — SARS-CoV-2 Wuhan-Hu-1 reference genome
#   references/NC_045512.2.gff3        — annotation, for SnpEff database build

set -euo pipefail

SAMPLE_ID="${1:?sample_id required}"
R1="${2:?trimmed R1 required}"
R2="${3:?trimmed R2 required}"
OUTDIR="${4:?outdir required}"
THREADS="${5:-8}"

REF="${REF_SARSCOV2:-references/NC_045512.2.fasta}"
MIN_COVERAGE=85       # percent
MIN_DEPTH=51          # X

mkdir -p "$OUTDIR/$SAMPLE_ID"
cd "$OUTDIR/$SAMPLE_ID"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "[02_virus_analysis] Indexing reference genome..."
    bwa index "$REF"
fi

echo "[02_virus_analysis] $SAMPLE_ID — BWA-MEM mapping to SARS-CoV-2 reference"
bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -o "${SAMPLE_ID}.sorted.bam" -
samtools index "${SAMPLE_ID}.sorted.bam"

echo "[02_virus_analysis] $SAMPLE_ID — Qualimap mapping QC"
qualimap bamqc -bam "${SAMPLE_ID}.sorted.bam" -outdir "qualimap_${SAMPLE_ID}" \
    --java-mem-size=4G

echo "[02_virus_analysis] $SAMPLE_ID — computing coverage/depth"
samtools depth -a "${SAMPLE_ID}.sorted.bam" > "${SAMPLE_ID}.depth.txt"
COVERAGE=$(awk '$3>0' "${SAMPLE_ID}.depth.txt" | wc -l)
GENOME_LEN=$(awk 'END{print NR}' "${SAMPLE_ID}.depth.txt")
PCT_COV=$(awk -v c="$COVERAGE" -v g="$GENOME_LEN" 'BEGIN{printf "%.2f", (c/g)*100}')
MEAN_DEPTH=$(awk '{sum+=$3} END{printf "%.1f", sum/NR}' "${SAMPLE_ID}.depth.txt")
echo "  Genome coverage: ${PCT_COV}%  |  Mean depth: ${MEAN_DEPTH}X"

echo "[02_virus_analysis] $SAMPLE_ID — de novo assembly with SPAdes"
spades.py --isolate \
    -1 "$R1" -2 "$R2" \
    -o "spades_${SAMPLE_ID}" \
    -t "$THREADS" -m 32

echo "[02_virus_analysis] $SAMPLE_ID — MAFFT alignment against reference"
cat "$REF" "spades_${SAMPLE_ID}/scaffolds.fasta" > "${SAMPLE_ID}_combined.fasta"
mafft --auto --thread "$THREADS" "${SAMPLE_ID}_combined.fasta" > "${SAMPLE_ID}_aligned.fasta"

echo "[02_virus_analysis] $SAMPLE_ID — Nextclade lineage/clade assignment"
nextclade run \
    --input-dataset references/nextclade_sars-cov-2_dataset \
    --output-all "nextclade_${SAMPLE_ID}" \
    "spades_${SAMPLE_ID}/scaffolds.fasta"

# Variant calling only proceeds if coverage/depth thresholds are met,
# matching the manuscript's reliability criteria for LoFreq calls.
PASS_COV=$(awk -v p="$PCT_COV" -v m="$MIN_COVERAGE" 'BEGIN{print (p>=m)?1:0}')
PASS_DEPTH=$(awk -v d="$MEAN_DEPTH" -v m="$MIN_DEPTH" 'BEGIN{print (d>m)?1:0}')

if [[ "$PASS_COV" -eq 1 && "$PASS_DEPTH" -eq 1 ]]; then
    echo "[02_virus_analysis] $SAMPLE_ID — QC thresholds met, calling variants with LoFreq"
    lofreq indelqual --dindel -f "$REF" -o "${SAMPLE_ID}.indelqual.bam" "${SAMPLE_ID}.sorted.bam"
    samtools index "${SAMPLE_ID}.indelqual.bam"
    lofreq call -f "$REF" -o "${SAMPLE_ID}.lofreq.vcf" "${SAMPLE_ID}.indelqual.bam"

    echo "[02_virus_analysis] $SAMPLE_ID — SnpEff annotation"
    snpEff -v NC_045512.2 "${SAMPLE_ID}.lofreq.vcf" > "${SAMPLE_ID}.annotated.vcf"
else
    echo "[02_virus_analysis] $SAMPLE_ID — SKIPPED variant calling (coverage ${PCT_COV}% / depth ${MEAN_DEPTH}X below thresholds ${MIN_COVERAGE}%/${MIN_DEPTH}X)"
fi

echo "[02_virus_analysis] $SAMPLE_ID — done"
