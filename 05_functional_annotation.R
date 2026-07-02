#!/usr/bin/env Rscript
#
# 05_functional_annotation.R — Functional annotation of DEGs (Figure 1, panel C, step 4)
#
# Queries the STRING database for protein-protein interaction (PPI) networks
# and the Enrichr API for Gene Ontology functional enrichment.
#
# Tools: STRING db v11.5 (min. interaction confidence score 0.4), Enrichr
#
# Usage:
#   Rscript 05_functional_annotation.R --deg_table <DEG_results.tsv> --outdir <dir>
#
# Requires: STRINGdb, enrichR (install via BiocManager::install("STRINGdb")
#           and install.packages("enrichR"))

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option("--deg_table", type = "character", help = "DEG_results.tsv from 04_deseq2_analysis.R"),
  make_option("--outdir", type = "character", default = "05_functional", help = "Output directory"),
  make_option("--species", type = "integer", default = 9606, help = "NCBI taxonomy ID [default 9606 = human]"),
  make_option("--string_score_threshold", type = "integer", default = 400,
              help = "STRING minimum combined score, 0-1000 scale [default 400 = confidence 0.4]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

deg <- read_tsv(opt$deg_table, show_col_types = FALSE)

if (nrow(deg) == 0) {
  stop("No DEGs found in ", opt$deg_table, " — nothing to annotate.")
}

gene_list <- unique(deg$gene_id)
write_lines(gene_list, file.path(opt$outdir, "deg_gene_list.txt"))

# --- STRING PPI network ---------------------------------------------------
if (requireNamespace("STRINGdb", quietly = TRUE)) {
  library(STRINGdb)

  string_db <- STRINGdb$new(
    version = "11.5",
    species = opt$species,
    score_threshold = opt$string_score_threshold,
    input_directory = ""
  )

  gene_df <- data.frame(gene = gene_list)
  mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

  write_tsv(mapped, file.path(opt$outdir, "string_mapped_genes.tsv"))

  png(file.path(opt$outdir, "string_ppi_network.png"), width = 1600, height = 1600, res = 150)
  string_db$plot_network(mapped$STRING_id)
  dev.off()

  # Hub genes by degree of connectivity
  interactions <- string_db$get_interactions(mapped$STRING_id)
  degree_tbl <- bind_rows(
    interactions %>% count(from, name = "n"),
    interactions %>% count(to, name = "n") %>% rename(from = to)
  ) %>%
    group_by(from) %>%
    summarise(degree = sum(n)) %>%
    arrange(desc(degree))
  write_tsv(degree_tbl, file.path(opt$outdir, "string_hub_genes.tsv"))

  cat("[05_functional_annotation] STRING network + hub genes written.\n")
} else {
  cat("[05_functional_annotation] STRINGdb package not installed — skipping PPI network.\n",
      "Install with: BiocManager::install('STRINGdb')\n")
}

# --- Enrichr functional enrichment ----------------------------------------
if (requireNamespace("enrichR", quietly = TRUE)) {
  library(enrichR)

  setEnrichrSite("Enrichr")
  dbs <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
           "KEGG_2021_Human", "Reactome_2022")

  enriched <- enrichr(gene_list, dbs)

  for (db_name in names(enriched)) {
    out_file <- file.path(opt$outdir, paste0("enrichr_", db_name, ".tsv"))
    write_tsv(enriched[[db_name]], out_file)
  }

  cat("[05_functional_annotation] Enrichr results written for:", paste(dbs, collapse = ", "), "\n")
} else {
  cat("[05_functional_annotation] enrichR package not installed — skipping enrichment.\n",
      "Install with: install.packages('enrichR')\n")
}

cat("[05_functional_annotation] Done.\n")
