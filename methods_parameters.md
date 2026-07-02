# Exact parameters and thresholds used in the manuscript

Reproduced here for transparency and reproducibility. These values are what the scripts in `scripts/` are configured to use by default.

## Sample cohort

- Initial extraction: 46 samples
- Final cohort after RNA-quality filtering (RNA yield, then RIN score within each group): **34 samples**
  - Severe: n = 9 (from 11)
  - Moderate: n = 14 (from 19)
  - Negative: n = 11 (from 16)
- RIN values, final cohort: Severe 1.4–6.2 (median 5.4); Moderate 1.0–5.2 (median 1.5)
- Sequencing: paired-end 100 bp, Illumina HiSeq X, ≥37.4 million read pairs/sample

## Preprocessing

- **FastQC** v0.12.0 — per-read QC
- **MultiQC** v1.14 — aggregated QC report
- **Cutadapt** v1.18 — adapter trimming (default Illumina adapter set; adjust to your library prep kit)

## Virus analysis

- Reference genome: SARS-CoV-2 isolate Wuhan-Hu-1, **NC_045512.2**
- Mapping: **BWA-MEM** v0.7.17 (Li & Durbin, 2009)
- Mapping QC: **Qualimap** v2.2.1
- Assembly: **SPAdes** v3.15.3
- Alignment for consensus/FASTA generation: **MAFFT** v7.471
- Lineage/clade assignment: **Nextclade** v1.5.3 (Nextstrain toolset)
- Variant calling: **LoFreq** v2.1.5
  - Minimum genome coverage: **85%**
  - Minimum depth: **>51X**
- Variant annotation: **SnpEff**
- Downstream stats/visualization: Python — pandas, matplotlib, seaborn

## Host transcriptional response

- Reference genome: **GRCh38**
- Mapping: **BWA-MEM** v0.7.17
- Quantification: **FeatureCounts** (Subread) v2.0.1
- Differential expression: **DESeq2** (R/Bioconductor)
  - Significance criteria (both required): |log2FoldChange| ≥ 1 (≥2-fold) **and** adjusted p-value (Benjamini–Hochberg) < 0.05
- PPI network: **STRING** database v11.5 (https://string-db.org)
  - Minimum interaction confidence score: **0.4**
- Functional annotation / pathway enrichment: **EnrichR** (https://maayanlab.cloud/Enrichr/)
  - Ranking: combined score integrating Fisher's exact test p-value and z-score

> **Note on rRNA depletion:** rRNA depletion removes structural rRNA (18S, 28S, 5S, mitochondrial rRNA) but does **not** deplete mRNAs encoding ribosomal proteins (RPS/RPL gene families). Differential expression detected for RPS/RPL transcripts therefore reflects genuine transcriptional changes, not a depletion artifact.

## Microbiome analysis

- Host/viral read removal: **Bowtie2** v2.4.5
- Chimera filtering: **VSEARCH** v2.22.0
- Taxonomic profiling: **MetaPhlAn** v4.0.6 (marker-gene database)
- Alpha diversity: **Vegan** package in R

## Ethics

- Institut Pasteur de Tunis Ethics Committee approval ID: **2021/20/I**
- Written informed consent obtained from all participants
- All data anonymized
- Negative extraction controls (VTM blanks) included in every extraction batch; confirmed RT-qPCR negative
