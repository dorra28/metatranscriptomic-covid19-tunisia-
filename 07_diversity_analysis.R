#!/usr/bin/env Rscript
#
# 07_diversity_analysis.R — Microbiome alpha diversity (Figure 1, panel D, step 4)
#
# Merges per-sample MetaPhlAn4 species-level abundance profiles into one
# abundance matrix and computes alpha diversity metrics with the Vegan package.
#
# Tools: Vegan (R/CRAN)
#
# Usage:
#   Rscript 07_diversity_analysis.R --abundance_dir <dir> --samples <samples.tsv> --outdir <dir>
#
# Expects one <sample_id>/<sample_id>_metaphlan_profile.tsv file per sample
# (produced by 06_microbiome_analysis.sh).

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(tidyverse)
})

option_list <- list(
  make_option("--abundance_dir", type = "character", help = "Directory containing per-sample MetaPhlAn output subfolders"),
  make_option("--samples", type = "character", help = "Sample sheet TSV (sample_id, group, ...)"),
  make_option("--outdir", type = "character", default = "07_diversity", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

samples <- read_tsv(opt$samples, show_col_types = FALSE)

# --- Read + merge per-sample MetaPhlAn species-level profiles -------------
read_metaphlan_species <- function(sample_id, abundance_dir) {
  f <- file.path(abundance_dir, sample_id, paste0(sample_id, "_metaphlan_profile.tsv"))
  if (!file.exists(f)) {
    warning("Missing MetaPhlAn profile for ", sample_id, " — skipping")
    return(NULL)
  }
  df <- read_tsv(f, comment = "#", col_names = c("clade_name", "clade_taxid", "relative_abundance", "additional_species"),
                  show_col_types = FALSE)
  df %>%
    filter(str_detect(clade_name, "s__")) %>%          # species-level rows only
    mutate(species = str_extract(clade_name, "s__[^|]+")) %>%
    select(species, relative_abundance) %>%
    rename(!!sample_id := relative_abundance)
}

profiles <- lapply(samples$sample_id, read_metaphlan_species, abundance_dir = opt$abundance_dir)
profiles <- profiles[!sapply(profiles, is.null)]

abundance_matrix <- reduce(profiles, full_join, by = "species") %>%
  column_to_rownames("species") %>%
  as.matrix()
abundance_matrix[is.na(abundance_matrix)] <- 0

write.csv(abundance_matrix, file.path(opt$outdir, "species_abundance_matrix.csv"))

# --- Alpha diversity (Vegan) ------------------------------------------------
# Samples as rows for vegan
mat_t <- t(abundance_matrix)

shannon <- diversity(mat_t, index = "shannon")
simpson <- diversity(mat_t, index = "simpson")
observed_richness <- specnumber(mat_t)
chao1 <- estimateR(round(mat_t))["S.chao1", ]

diversity_df <- tibble(
  sample_id = rownames(mat_t),
  shannon = shannon,
  simpson = simpson,
  observed_richness = observed_richness,
  chao1 = chao1
) %>%
  left_join(samples %>% select(sample_id, group), by = "sample_id")

write_tsv(diversity_df, file.path(opt$outdir, "alpha_diversity.tsv"))

# --- Group comparison (Kruskal-Wallis, non-parametric across 3 groups) ----
metrics <- c("shannon", "simpson", "observed_richness", "chao1")
kw_results <- map_dfr(metrics, function(m) {
  test <- kruskal.test(diversity_df[[m]] ~ diversity_df$group)
  tibble(metric = m, kruskal_statistic = test$statistic, p_value = test$p.value)
})
write_tsv(kw_results, file.path(opt$outdir, "alpha_diversity_group_comparison.tsv"))

# --- Boxplots ----------------------------------------------------------
pdf(file.path(opt$outdir, "alpha_diversity_boxplots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (m in metrics) {
  boxplot(diversity_df[[m]] ~ diversity_df$group,
          main = m, xlab = "Group", ylab = m, col = c("#4C72B0", "#DD8452", "#55A868"))
}
dev.off()

cat("[07_diversity_analysis] Done. Alpha diversity metrics for", nrow(diversity_df), "samples written to",
    file.path(opt$outdir, "alpha_diversity.tsv"), "\n")
