# Uncovering COVID-19 Dynamics in Tunisian Patients: A Meta-Transcriptomic Approach

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21128187.svg)](https://doi.org/10.5281/zenodo.21128187)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible bioinformatics pipeline accompanying the manuscript *"Uncovering COVID-19 Dynamics in Tunisian Patients: A Meta-Transcriptomic Approach"* (submitted to *Scientific Reports*).

> **Authors:** Dorra Rjaibi, Oussama Souiai, Lilia Romdhane — Institut Pasteur de Tunis

This repository contains the scripts used to process paired-end metatranscriptomic sequencing data from nasopharyngeal swabs of Tunisian COVID-19 patients (severe, moderate, and RT-PCR-negative controls) for three parallel analyses:

1. **Virus analysis** — SARS-CoV-2 genome mapping, assembly, lineage/clade assignment, and variant calling
2. **Host transcriptional response** — differential gene expression and functional/pathway enrichment
3. **Microbiome profiling** — taxonomic classification and alpha diversity of the nasopharyngeal microbiome

---

## Study overview

- **Samples:** 34 Tunisian nasopharyngeal swab (NPS) samples (Severe n=9, Moderate n=14, Negative n=11), collected November 2021–January 2022 (Omicron-dominant period)
- **Ethics approval:** Institut Pasteur de Tunis Ethics Committee, ID 2021/20/I
- **RNA extraction:** QIAamp Viral RNA Mini Kit (Qiagen), 140 µL input
- **Sequencing:** Paired-end 100 bp, Illumina HiSeq X, ≥37.4 million read pairs/sample

## Pipeline overview

```
                         ┌────────────────────┐
                         │   Raw FASTQ files   │
                         └─────────┬───────────┘
                                   │
                    ┌──────────────▼───────────────┐
                    │  A. Preprocessing              │
                    │  FastQC → MultiQC → Cutadapt   │
                    └──────────────┬───────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        │                          │                           │
┌───────▼────────┐      ┌──────────▼──────────┐     ┌──────────▼──────────┐
│ B. Virus        │      │ C. Host response     │     │ D. Microbiome        │
│ analysis        │      │                       │     │ profiling             │
│                  │      │                       │     │                       │
│ BWA-MEM →        │      │ BWA-MEM (GRCh38) →   │     │ Bowtie2 (host/viral   │
│ Samtools →       │      │ Samtools →            │     │ removal) → VSEARCH    │
│ Qualimap →       │      │ FeatureCounts →       │     │ (chimera filtering) → │
│ SPAdes →         │      │ DESeq2 →              │     │ MetaPhlAn4 →          │
│ MAFFT →          │      │ EnrichR + STRING      │     │ Vegan (alpha          │
│ Nextclade →      │      │                       │     │ diversity, R)         │
│ LoFreq → SnpEff  │      │                       │     │                       │
└──────────────────┘      └───────────────────────┘     └───────────────────────┘
```

This matches **Figure 1** of the manuscript (Optimized metatranscriptomic pipeline for SARS-CoV-2, host and microbiome analyses).

## Repository structure

```
.
├── README.md
├── LICENSE
├── environment.yml            # conda environment with pinned tool versions
├── config/
│   └── samples_template.tsv   # sample sheet template (fill in with your sample IDs / paths)
├── scripts/
│   ├── 00_run_pipeline.sh     # master script, calls everything below in order
│   ├── 01_qc_trimming.sh      # FastQC + MultiQC + Cutadapt
│   ├── 02_virus_analysis.sh   # BWA-MEM, Qualimap, SPAdes, MAFFT, Nextclade, LoFreq, SnpEff
│   ├── 03_host_mapping_quant.sh # BWA-MEM (GRCh38) + FeatureCounts
│   ├── 04_deseq2_analysis.R   # DESeq2 differential expression
│   ├── 05_functional_annotation.R # STRING + EnrichR enrichment
│   ├── 06_microbiome_analysis.sh  # Bowtie2, VSEARCH, MetaPhlAn4
│   └── 07_diversity_analysis.R    # Vegan alpha diversity
└── docs/
    └── methods_parameters.md  # exact tool versions & thresholds used in the manuscript
```

## Tool versions (as reported in Methods)

| Step | Tool | Version |
|---|---|---|
| QC | FastQC | 0.12.0 |
| QC aggregation | MultiQC | 1.14 |
| Adapter trimming | Cutadapt | 1.18 |
| Read mapping (virus & host) | BWA-MEM | 0.7.17 |
| Mapping QC | Qualimap | 2.2.1 |
| Assembly | SPAdes | 3.15.3 |
| Alignment | MAFFT | 7.471 |
| Lineage/clade assignment | Nextclade | 1.5.3 |
| Variant calling | LoFreq | 2.1.5 |
| Variant annotation | SnpEff | latest |
| Gene quantification | FeatureCounts (Subread) | 2.0.1 |
| Differential expression | DESeq2 (R/Bioconductor) | Bioconductor release |
| PPI network | STRING db | v11.5 |
| Functional enrichment | EnrichR | web API |
| Host/viral read removal | Bowtie2 | 2.4.5 |
| Chimera filtering | VSEARCH | 2.22.0 |
| Taxonomic profiling | MetaPhlAn | 4.0.6 |
| Alpha diversity | Vegan (R) | CRAN release |

See `docs/methods_parameters.md` for the exact thresholds (coverage/depth cutoffs, log2FC/padj cutoffs, STRING confidence score, etc.) reproduced from the manuscript.

## Requirements

- Linux (tested on Ubuntu 20.04+)
- [conda](https://docs.conda.io/) / [mamba](https://mamba.readthedocs.io/) for environment management
- R ≥ 4.2 with Bioconductor ≥ 3.16
- ~16 GB RAM minimum recommended for SPAdes assembly and MetaPhlAn database loading

## Quick start

```bash
git clone https://github.com/<your-username>/metatranscriptomic-covid19-tunisia.git
cd metatranscriptomic-covid19-tunisia

# create the environment
conda env create -f environment.yml
conda activate covid19-metatx

# fill in config/samples_template.tsv with your sample IDs and FASTQ paths
cp config/samples_template.tsv config/samples.tsv
nano config/samples.tsv

# run the full pipeline
bash scripts/00_run_pipeline.sh config/samples.tsv
```

Each script can also be run independently — see the header comments in each file for required inputs/outputs.

## Citation

If you use this pipeline, please cite both the software and the manuscript:

**Software:**
> Rjaibi D, Souiai O., Romdhane L (2026). *Metatranscriptomic COVID-19 Tunisia: Pipeline for SARS-CoV-2, Host, and Microbiome Analysis* (v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.21128187

**Manuscript:**
> Rjaibi D et al .2026  Uncovering COVID-19 Dynamics in Tunisian Patients: A Meta-Transcriptomic Approach. *Scientific Reports* (submitted).


## Ethics statement

This study was approved by the Ethics Committee of the Institut Pasteur de Tunis (ID 2021/20/I). Informed consent was obtained from all participants; all data were anonymized.

## License

MIT License — see [LICENSE](LICENSE).
