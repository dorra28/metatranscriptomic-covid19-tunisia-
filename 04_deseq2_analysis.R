#!/usr/bin/env Rscript
#
# 04_deseq2_analysis.R — Host differential gene expression (Figure 1, panel C, step 3)
#
# Merges per-sample FeatureCounts output, runs DESeq2, and retains genes
# passing the manuscript's significance criteria:
#   |log2FoldChange| >= 1 (>= 2-fold change) AND
#   adjusted p-value (Benjamini-Hochberg) < 0.05
#
# Usage:
#   Rscript 04_deseq2_analysis.R --counts_dir <dir> --samples <samples.tsv> --outdir <dir>
#
# Expects one <sample_id>.counts.tsv file per sample in counts_dir
# (produced by 03_host_mapping_quant.sh), with columns: gene_id, count.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(tidyverse)
})

option_list <- list(
  make_option("--counts_dir", type = "character", help = "Directory with per-sample *.counts.tsv files"),
  make_option("--samples", type = "character", help = "Sample sheet TSV (sample_id, group, ...)"),
  make_option("--outdir", type = "character", default = "04_deg", help = "Output directory"),
  make_option("--lfc_threshold", type = "double", default = 1, help = "|log2FC| threshold [default 1]"),
  make_option("--padj_threshold", type = "double", default = 0.05, help = "Adjusted p-value threshold [default 0.05]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Load sample sheet -------------------------------------------------
samples <- read_tsv(opt$samples, show_col_types = FALSE)

# --- Merge per-sample count files into one count matrix -----------------
count_files <- file.path(opt$counts_dir, paste0(samples$sample_id, ".counts.tsv"))
missing <- count_files[!file.exists(count_files)]
if (length(missing) > 0) {
  stop("Missing count files: ", paste(missing, collapse = ", "))
}

count_list <- lapply(seq_along(count_files), function(i) {
  df <- read_tsv(count_files[i], col_names = c("gene_id", samples$sample_id[i]), show_col_types = FALSE)
  df
})

count_matrix <- reduce(count_list, full_join, by = "gene_id") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
count_matrix[is.na(count_matrix)] <- 0
storage.mode(count_matrix) <- "integer"

# --- Build coldata / metadata -------------------------------------------
coldata <- samples %>%
  column_to_rownames("sample_id")
coldata <- coldata[colnames(count_matrix), , drop = FALSE]
coldata$group <- factor(coldata$group, levels = c("Negative", "Moderate", "Severe"))

# --- DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ group
)

# Filter very low-count genes prior to testing
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Pairwise contrasts vs Negative controls (adjust as needed for your design)
contrasts <- list(
  Severe_vs_Negative   = c("group", "Severe", "Negative"),
  Moderate_vs_Negative = c("group", "Moderate", "Negative"),
  Severe_vs_Moderate   = c("group", "Severe", "Moderate")
)

all_results <- list()

for (contrast_name in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[contrast_name]], alpha = opt$padj_threshold)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(contrast = contrast_name)

  write_tsv(res_df, file.path(opt$outdir, paste0(contrast_name, "_full_results.tsv")))

  sig_df <- res_df %>%
    filter(!is.na(padj), padj < opt$padj_threshold, abs(log2FoldChange) >= opt$lfc_threshold)

  write_tsv(sig_df, file.path(opt$outdir, paste0(contrast_name, "_DEG_significant.tsv")))

  all_results[[contrast_name]] <- sig_df
}

combined_deg <- bind_rows(all_results)
write_tsv(combined_deg, file.path(opt$outdir, "DEG_results.tsv"))

# Normalized counts, for downstream visualization
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")
write_tsv(norm_counts, file.path(opt$outdir, "normalized_counts.tsv"))

cat(sprintf(
  "\n[04_deseq2_analysis] Done. %d significant DEGs across %d contrasts (|log2FC|>=%.1f, padj<%.2f)\n",
  nrow(combined_deg), length(contrasts), opt$lfc_threshold, opt$padj_threshold
))
