#!/usr/bin/env Rscript
# Shared data loading + helpers for poster figures (CNPN 2026).
# Sources the main loader, derives sex, builds batch-level matrices,
# defines theme/palettes, and provides a permutation-validated PLS-DA helper.

# --- Source the standard loader from the parent dir ------------------------
.orig_wd <- getwd()
setwd("..")
source("01_database_analysis.R")
setwd(.orig_wd)

suppressPackageStartupMessages({
  library(cowplot)
  library(patchwork)
  library(ggrepel)
  library(viridis)
})

# --- Sex derivation ---------------------------------------------------------
peptides <- peptides %>%
  mutate(
    sex = case_when(
      grepl("_female(\\b|$)", batch, ignore.case = TRUE) ~ "Female",
      grepl("_male(\\b|$)",   batch, ignore.case = TRUE) ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female"))
  )

cat("\n--- Sex parsing ---\n")
print(peptides %>% count(batch, tissue, sex))

# --- Palettes ---------------------------------------------------------------
SEX_COLORS    <- c(Male = "#1f77b4", Female = "#e377c2")
TISSUE_COLORS <- c(
  Cerebrum       = "#1f78b4",
  BrainStem      = "#33a02c",
  Cerebellum     = "#e31a1c",
  MockWholeBrain = "#6a3d9a"
)
DIRECTION_COLORS <- c(`Up in Male` = "#1f77b4", `Up in Female` = "#e377c2", `n.s.` = "grey80")

# --- Phosphosite parser ----------------------------------------------------
# Extracts site notation (e.g. "S14", "T8,S22") from a modified peptide
# sequence. Position is RELATIVE TO THE PEPTIDE (not the protein) because
# protein-coordinate sites need a UniProt sequence lookup we don't have.
# Examples:
#   "LSEGPAALAGPAS[79.9663]PPR"      -> "S14"
#   "TEEVLSPDGS[79.9663]PSKS[79.9663]PSK" -> "S10,S14"
#   "S[79.9663]PGM[15.9949]LEPLGSAR" -> "S1"
parse_phosphosite <- function(modseq) {
  if (is.na(modseq) || modseq == "") return(NA_character_)
  # Locate every "[79.9663]" marker in the modified string
  matches <- gregexpr("\\[79\\.9663\\]", modseq)[[1]]
  if (matches[1] == -1) return(NA_character_)
  ml <- attr(matches, "match.length")
  sites <- vapply(seq_along(matches), function(i) {
    pre  <- substr(modseq, 1, matches[i] - 1)
    bare <- gsub("\\[[^\\]]*\\]", "", pre)  # strip ALL prior modifications
    res  <- substr(bare, nchar(bare), nchar(bare))
    paste0(res, nchar(bare))
  }, character(1))
  paste(sites, collapse = ",")
}

# Vectorized helper: gene + " p" + sites — e.g. "Lrrc7 pS14"
phosphosite_label <- function(gene, modseq) {
  ifelse(is.na(gene) | gene == "", NA_character_,
         vapply(seq_along(modseq), function(i) {
           sites <- parse_phosphosite(modseq[i])
           if (is.na(sites)) gene[i]
           else sprintf("%s p%s", gene[i], sites)
         }, character(1)))
}

# --- QC drop-out flag -------------------------------------------------------
# Cerebellum_TCPG1_female had ~7,899 rows (vs ~140k for sibling batches): clear
# upstream QC failure. Excluded from multivariate analyses (PCA/PLS-DA).
MV_BAD_BATCHES <- c("Cerebellum_TCPG1_female")

# --- Save destination -------------------------------------------------------
# Poster-ready figures: clean folder containing only FIG1-FIG4 + FIG3b.
# Legacy outputs (P1a/F0x/older composites) remain in ../figures_poster.
FIG_DIR <- "../sysq_mouse_figs"
TAB_DIR <- "../tables"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Theme ------------------------------------------------------------------
theme_poster <- function(base_size = 13) {
  theme_publication(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 3, face = "bold"),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30"),
      plot.caption  = element_text(size = base_size - 3, color = "grey45", hjust = 0),
      axis.title    = element_text(size = base_size,     face = "bold"),
      axis.text     = element_text(size = base_size - 2),
      legend.title  = element_text(size = base_size - 1, face = "bold"),
      legend.text   = element_text(size = base_size - 2),
      strip.text    = element_text(size = base_size - 1, face = "bold"),
      panel.grid.major = element_line(size = 0.25, color = "grey92"),
      panel.grid.minor = element_blank()
    )
}

save_poster <- function(plot, name, width = 14, height = 10) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}

# --- Batch-median normalization on log10(heavy_intensity) ------------------
# Female batches have ~2x higher raw heavy intensity than male siblings — a
# loading/MS effect, not biology. Apply the same batch-median normalization
# that the human SysQuan paper uses (subtract per-batch median of
# log10(heavy), add back the global median). The normalized column is
# carried as `log_heavy_norm` and is the canonical signal for downstream DE.
.global_log_heavy <- peptides %>%
  filter(!is.na(heavy_intensity), heavy_intensity > 0) %>%
  pull(heavy_intensity) %>% log10() %>% median(na.rm = TRUE)

.batch_log_heavy_med <- peptides %>%
  filter(!is.na(heavy_intensity), heavy_intensity > 0) %>%
  group_by(batch) %>%
  summarise(batch_med = median(log10(heavy_intensity), na.rm = TRUE),
            .groups = "drop")

cat("\n--- Batch-median normalization (log10 heavy) ---\n")
cat(sprintf("Global median log10(heavy) = %.3f\n", .global_log_heavy))
print(.batch_log_heavy_med)

peptides <- peptides %>%
  left_join(.batch_log_heavy_med, by = "batch") %>%
  mutate(
    log10_heavy     = ifelse(!is.na(heavy_intensity) & heavy_intensity > 0,
                             log10(heavy_intensity), NA_real_),
    log_heavy_norm  = log10_heavy - batch_med + .global_log_heavy
  ) %>%
  select(-batch_med)

# --- Batch metadata ---------------------------------------------------------
batch_meta <- peptides %>%
  select(batch, tissue, sex) %>%
  distinct() %>%
  filter(!is.na(sex)) %>%
  mutate(group = paste(tissue, sex, sep = " · "))

# --- Batch-level protein intensity matrix ----------------------------------
# Mean log_heavy_norm (batch-median-normalized log10 heavy) per protein
# per batch. Used by PCA / PLS-DA.
build_protein_matrix <- function(df) {
  df %>%
    filter(!is.na(log_heavy_norm), !is.na(sex)) %>%
    group_by(protein_id, batch) %>%
    summarise(mean_log10 = mean(log_heavy_norm, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = protein_id, values_from = mean_log10, values_fill = NA)
}

batch_matrix_full <- build_protein_matrix(peptides)
batch_matrix_phos <- build_protein_matrix(peptides %>% filter(phospho))

cat(sprintf("Full protein matrix: %d batches x %d proteins\n",
            nrow(batch_matrix_full), ncol(batch_matrix_full) - 1))
cat(sprintf("Phospho matrix:      %d batches x %d proteins\n",
            nrow(batch_matrix_phos), ncol(batch_matrix_phos) - 1))

# --- Helpers: matrix preparation for mixOmics PLS-DA -----------------------
# `meta`        a tibble with columns `batch` and `group` (factor)
# Returns a list with X (samples × features, NA→feature-mean), Y aligned to
# rownames(X), and the meta subset corresponding to the kept batches.
prep_plsda <- function(batch_mat, meta, drop_batches = MV_BAD_BATCHES,
                       top_var = 2000) {
  stopifnot(all(c("batch", "group") %in% colnames(meta)))
  X <- batch_mat %>%
    dplyr::filter(!batch %in% drop_batches) %>%
    tibble::column_to_rownames("batch") %>%
    as.matrix()
  # NA → column (feature) mean
  for (j in seq_len(ncol(X))) {
    m <- mean(X[, j], na.rm = TRUE)
    X[is.na(X[, j]), j] <- m
  }
  # Drop zero-variance features
  vars <- apply(X, 2, var, na.rm = TRUE)
  X <- X[, !is.na(vars) & vars > 0, drop = FALSE]
  # Top-variance feature filter
  if (!is.null(top_var) && ncol(X) > top_var) {
    vars <- apply(X, 2, var)
    X <- X[, order(-vars)[seq_len(top_var)], drop = FALSE]
  }
  # Align Y to rownames(X) by batch ID
  meta_aligned <- meta[match(rownames(X), meta$batch), , drop = FALSE]
  Y <- droplevels(factor(meta_aligned$group))
  list(X = X, Y = Y, batches = rownames(X), meta = meta_aligned)
}

# --- Permutation-validated PLS-DA ------------------------------------------
# Returns: ncomp=2 sPLS-DA scores, permutation p-value (label shuffle, n_perm
# times, comparing observed between-group / within-group SS ratio to null).
plsda_with_perm <- function(X, Y, n_perm = 999, ncomp = 2, seed = 1) {
  stopifnot(requireNamespace("mixOmics", quietly = TRUE))
  set.seed(seed)
  fit <- mixOmics::plsda(X, Y, ncomp = ncomp, scale = TRUE)
  scores <- as.data.frame(fit$variates$X)
  colnames(scores) <- paste0("LV", seq_len(ncol(scores)))
  scores$batch <- rownames(scores)
  scores$Y <- Y

  # Test statistic: between/within SS on LV1+LV2
  ss <- function(scores_mat, labels) {
    centers <- t(sapply(split(seq_len(nrow(scores_mat)), labels),
                        function(idx) colMeans(scores_mat[idx, , drop = FALSE])))
    grand <- colMeans(scores_mat)
    bw <- sum(table(labels) * rowSums((centers - matrix(grand, nrow = nrow(centers),
                                                        ncol = length(grand), byrow = TRUE))^2))
    wi <- 0
    for (lab in unique(labels)) {
      idx <- which(labels == lab)
      wi <- wi + sum((scores_mat[idx, , drop = FALSE] -
                      matrix(centers[as.character(lab), ], nrow = length(idx),
                             ncol = ncol(scores_mat), byrow = TRUE))^2)
    }
    bw / max(wi, 1e-9)
  }
  obs <- ss(as.matrix(scores[, paste0("LV", seq_len(ncomp))]), Y)

  # Permutation null
  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    Yp <- sample(Y)
    fp <- mixOmics::plsda(X, Yp, ncomp = ncomp, scale = TRUE)
    sp <- as.matrix(as.data.frame(fp$variates$X)[, seq_len(ncomp)])
    null[i] <- ss(sp, Yp)
  }
  p_val <- (sum(null >= obs) + 1) / (n_perm + 1)

  # Variance explained by each LV (mixOmics 'prop_expl_var')
  ve <- fit$prop_expl_var$X
  list(fit = fit, scores = scores, p_val = p_val, n_perm = n_perm,
       obs_stat = obs, null = null, var_explained = ve)
}

# --- Helpers: peptide-level Welch t-test for sex contrast within tissue ---
# `value_col` defaults to the batch-median-normalized heavy intensity.
# `channel` is kept for backward compat: if the named column doesn't already
# exist on `df`, we fall back to log10(channel) computed on the fly.
sex_de_within_tissue <- function(df, tissue_name,
                                 value_col = "log_heavy_norm",
                                 channel = "heavy_intensity",
                                 min_n_per_sex = 3) {
  if (value_col %in% colnames(df)) {
    d <- df %>%
      filter(tissue == tissue_name, !is.na(sex), !is.na(.data[[value_col]])) %>%
      mutate(log_int = .data[[value_col]])
  } else {
    d <- df %>%
      filter(tissue == tissue_name, !is.na(sex), !is.na(.data[[channel]]),
             .data[[channel]] > 0) %>%
      mutate(log_int = log10(.data[[channel]]))
  }
  pep_lvl <- d %>%
    group_by(protein_id, gene, sex, modified_peptide_sequence) %>%
    summarise(value = median(log_int, na.rm = TRUE), .groups = "drop")

  pep_lvl %>%
    group_by(protein_id, gene) %>%
    summarise(
      n_M = sum(sex == "Male"),
      n_F = sum(sex == "Female"),
      mean_M = mean(value[sex == "Male"], na.rm = TRUE),
      mean_F = mean(value[sex == "Female"], na.rm = TRUE),
      log2FC = (mean_M - mean_F) * log2(10),  # log10 → log2
      p_value = tryCatch({
        if (n_M >= min_n_per_sex && n_F >= min_n_per_sex)
          t.test(value[sex == "Male"], value[sex == "Female"], var.equal = FALSE)$p.value
        else NA_real_
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      adj_p = p.adjust(p_value, method = "BH"),
      tissue = tissue_name
    )
}

# --- GO BP enrichment from a gene-symbol vector vs background universe ----
# Lazy-loads clusterProfiler + org.Mm.eg.db on first use.
.universe_cache <- NULL
get_universe_syms <- function() {
  if (is.null(.universe_cache)) {
    .universe_cache <<- peptides %>%
      filter(!is.na(gene), gene != "") %>%
      pull(gene) %>% unique()
  }
  .universe_cache
}

enrich_gene_set <- function(gene_syms, label, fdr_cap = 0.10, top_n = 6,
                            universe = NULL, ont = "BP") {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })
  gene_syms <- unique(gene_syms[!is.na(gene_syms) & gene_syms != ""])
  if (length(gene_syms) < 3) {
    return(tibble(label = label, ontology = ont, Description = NA_character_,
                  p.adjust = NA_real_, gene_ratio_n = NA_real_,
                  Count = NA_integer_))
  }
  if (is.null(universe)) universe <- get_universe_syms()
  ego <- tryCatch(
    enrichGO(gene = gene_syms, universe = universe,
             OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = fdr_cap,
             readable = FALSE, minGSSize = 5, maxGSSize = 1000),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(tibble(label = label, ontology = ont, Description = NA_character_,
                  p.adjust = NA_real_, gene_ratio_n = NA_real_,
                  Count = NA_integer_))
  }
  as.data.frame(ego) %>%
    mutate(label = label, ontology = ont,
           gene_ratio_n = as.numeric(sub("/.*", "", GeneRatio)) /
                          as.numeric(sub(".*/", "", GeneRatio))) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    dplyr::select(label, ontology, Description, p.adjust, gene_ratio_n, Count)
}

# Run BP + MF + CC and return a single df with an ontology column
enrich_gene_set_all <- function(gene_syms, label, fdr_cap = 0.25, top_n = 5,
                                universe = NULL) {
  bind_rows(lapply(c("BP", "MF", "CC"), function(o) {
    enrich_gene_set(gene_syms, label, fdr_cap = fdr_cap, top_n = top_n,
                    universe = universe, ont = o)
  }))
}

# Standard 3-ontology dotplot used in FIG3b/FIG3c/FIG4. Faceted by ontology.
# Uses the original (unabbreviated) GO descriptions and wraps at the width
# that fits the panel cleanly — no in-place shortening (which made labels
# noisier than the wrap they were trying to avoid).
plot_go_3panels <- function(go_df, title, subtitle = NULL,
                            x_levels = c("Up in Male", "Up in Female"),
                            wrap_width = 42) {
  go_df <- go_df %>% filter(!is.na(Description))
  if (nrow(go_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 5, color = "grey40",
                      label = "No GO terms passed threshold\n(input lists too small)") +
             theme_void(base_size = 14) +
             labs(title = title, subtitle = subtitle) +
             theme(plot.background = element_rect(fill = "white", colour = NA),
                   plot.title    = element_text(face = "bold", size = 14,
                                                hjust = 0.5),
                   plot.subtitle = element_text(size = 11, color = "grey30",
                                                hjust = 0.5)))
  }
  go_df <- go_df %>%
    mutate(ontology = factor(ontology, levels = c("BP", "MF", "CC"),
                             labels = c("Biological Process",
                                        "Molecular Function",
                                        "Cellular Component")))
  ggplot(go_df,
         aes(x = factor(label, levels = x_levels),
             y = stringr::str_wrap(Description, wrap_width),
             size = Count, color = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ ontology, ncol = 1, scales = "free_y", strip.position = "top") +
    scale_color_gradientn(colors = GO_COLORS,
                          name = expression(-log[10]~"adj.p")) +
    scale_size_continuous(range = c(2.5, 8), name = "Gene\noverlap") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_poster() +
    theme(axis.text.y      = element_text(size = 10, lineheight = 0.95,
                                          color = "grey15"),
          axis.text.x      = element_text(size = 12),
          strip.background = element_rect(fill = "grey95", color = "grey70"),
          strip.text       = element_text(face = "bold", size = 12,
                                          margin = margin(3, 0, 3, 0)),
          panel.spacing.y  = unit(10, "pt"),
          legend.position  = "right",
          legend.box       = "vertical",
          plot.subtitle    = element_text(size = 10, color = "grey30"))
}

# Standard 3-color sequential color scale used across all GO dotplots
GO_COLORS <- c("#1f77b4", "#2ca02c", "#d62728")  # blue → green → red

# --- One-vs-rest tissue marker test ----------------------------------------
tissue_de_one_vs_rest <- function(df, target_tissue,
                                  value_col = "log_heavy_norm",
                                  channel   = "heavy_intensity",
                                  min_n = 3) {
  if (value_col %in% colnames(df)) {
    d <- df %>%
      filter(!is.na(.data[[value_col]])) %>%
      mutate(log_int = .data[[value_col]],
             grp = ifelse(tissue == target_tissue, "in", "out"))
  } else {
    d <- df %>%
      filter(!is.na(.data[[channel]]), .data[[channel]] > 0) %>%
      mutate(log_int = log10(.data[[channel]]),
             grp = ifelse(tissue == target_tissue, "in", "out"))
  }
  pep_lvl <- d %>%
    group_by(protein_id, gene, grp, modified_peptide_sequence) %>%
    summarise(value = median(log_int, na.rm = TRUE), .groups = "drop")

  pep_lvl %>%
    group_by(protein_id, gene) %>%
    summarise(
      n_in  = sum(grp == "in"),
      n_out = sum(grp == "out"),
      mean_in  = mean(value[grp == "in"],  na.rm = TRUE),
      mean_out = mean(value[grp == "out"], na.rm = TRUE),
      log2FC = (mean_in - mean_out) * log2(10),
      p_value = tryCatch({
        if (n_in >= min_n && n_out >= min_n)
          t.test(value[grp == "in"], value[grp == "out"], var.equal = FALSE)$p.value
        else NA_real_
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(adj_p = p.adjust(p_value, method = "BH"),
           tissue = target_tissue)
}
