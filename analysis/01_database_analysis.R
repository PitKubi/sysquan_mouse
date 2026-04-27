#!/usr/bin/env Rscript
# SysQuan Mouse Database Analysis - Data loading + preprocessing
# Filters the all-species bucket cache to mouse rows only.
# Generated 2026-04-23 from the human analysis in sysquant_export/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(scales)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(data.table)
})

# --- Publication theme ------------------------------------------------------
theme_publication <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(size = 0.3, color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(size = 1, color = "black"),
      axis.line        = element_line(size = 0.5, color = "black"),
      axis.ticks       = element_line(size = 0.5, color = "black"),
      axis.text        = element_text(size = 10, color = "black"),
      axis.title       = element_text(size = 12, face = "bold"),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.position  = "right",
      strip.background = element_rect(fill = "grey95", color = "black", size = 1),
      strip.text       = element_text(face = "bold", size = 10)
    )
}
theme_set(theme_publication())

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# --- Mouse tissue definition -----------------------------------------------
# Derived from production DB (species='mouse'): 4 brain regions, no overlap
# with human tissue labels.
MOUSE_TISSUES <- c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain")

# Stable palette for the four mouse tissues (used across figures)
TISSUE_COLORS <- c(
  Cerebrum       = "#1f78b4",
  BrainStem      = "#33a02c",
  Cerebellum     = "#e31a1c",
  MockWholeBrain = "#6a3d9a"
)

cat("=== SysQuan Mouse Database Analysis ===\n")

# --- Load bucket cache (all species), keep mouse only -----------------------
bucket_files <- list.files("../", pattern = "peptide_bucket_cache.*\\.csv$", full.names = TRUE)
# Prefer the all-species export (it contains the Phospho column)
has_phospho <- function(path) {
  header <- readLines(path, n = 1)
  grepl("Phospho", header, fixed = TRUE)
}
keep <- vapply(bucket_files, has_phospho, logical(1))
bucket_files <- bucket_files[keep]
if (!length(bucket_files)) stop("No bucket cache CSV with Phospho column found in ../")
latest_bucket <- bucket_files[order(file.mtime(bucket_files), decreasing = TRUE)][1]
cat("Reading:", latest_bucket, "\n")

# fread is ~10x faster than read_csv for a 400 MB file
peptides_dt <- fread(latest_bucket, showProgress = FALSE)

# Column names from the export command (kept aligned with the header order)
setnames(peptides_dt, c(
  "batch", "tissue", "instrument", "protein_id", "gene",
  "modified_peptide_sequence", "clean_peptide",
  "incorporation_pct", "doublet", "methionine", "cystein", "phospho",
  "proteotypic", "missed_cleavages", "length_valid",
  "hydrophobicity", "light_intensity", "heavy_intensity",
  "mrm_available", "prm_pasef_available", "mouse_std_quant",
  "median_mouse_conc", "mouse_lloq", "mouse_std_cv",
  "synthetic_peptide", "approved", "created_at"
))

# Filter to mouse tissues. Species is not stored in the CSV, but the four
# mouse tissue labels above are unique to mouse batches in production.
peptides_dt <- peptides_dt[tissue %in% MOUSE_TISSUES]

# Type cleanup
peptides_dt[, `:=`(
  light_intensity    = as.numeric(light_intensity),
  heavy_intensity    = as.numeric(heavy_intensity),
  hydrophobicity     = as.numeric(hydrophobicity),
  incorporation_pct  = suppressWarnings(as.numeric(incorporation_pct)),
  doublet            = as.logical(doublet),
  methionine         = as.logical(methionine),
  cystein            = as.logical(cystein),
  phospho            = as.logical(phospho),
  proteotypic        = as.logical(proteotypic),
  length_valid       = as.logical(length_valid),
  approved           = as.logical(approved)
)]

# Fix tissue factor order for consistent plotting
peptides_dt[, tissue := factor(tissue, levels = MOUSE_TISSUES)]

# Convert to tibble for tidyverse downstream; keep a DT alias in case we
# need fast operations later.
peptides <- as_tibble(peptides_dt)

cat("\nMouse data loaded!\n")
cat("  Rows:      ", nrow(peptides), "\n")
cat("  Proteins:  ", n_distinct(peptides$protein_id), "\n")
cat("  Batches:   ", n_distinct(peptides$batch), "\n")
cat("  Tissues:   ", paste(levels(peptides$tissue), collapse = ", "), "\n")
cat("  Phospho:   ", sum(peptides$phospho, na.rm = TRUE), "peptides\n")
cat("  Doublets:  ", sum(peptides$doublet, na.rm = TRUE), "\n\n")

tissue_summary <- peptides %>%
  group_by(tissue) %>%
  summarise(
    n_peptides      = n(),
    n_proteins      = n_distinct(protein_id),
    n_proteotypic   = sum(proteotypic, na.rm = TRUE),
    n_doublets      = sum(doublet, na.rm = TRUE),
    n_phospho       = sum(phospho, na.rm = TRUE),
    median_light_int = median(light_intensity, na.rm = TRUE),
    median_heavy_int = median(heavy_intensity, na.rm = TRUE),
    .groups = "drop"
  )
print(tissue_summary)
write_csv(tissue_summary, "tables/tissue_summary.csv")

cat("\n=== Generating Figures ===\n\n")
