#!/usr/bin/env Rscript
# FIGURE 2 — Phospho layer (Cerebrum, BrainStem, Cerebellum) (CNPN 2026)
# Panels:
#   A. PCA of phospho-protein intensity (Cer M/F, BS M/F, Cbl M; Cbl-F dropped)
#   B. PLS-DA on phospho-protein matrix with 999-perm validation p-value
#   C. Top phospho-protein region markers (heatmap, one-vs-rest, z-scored)
#   D. Bulk vs phospho coverage per region (Venn-style)
#
# Run from analysis/poster/:  Rscript FIG2_phospho.R

source("00_common.R")
suppressPackageStartupMessages({
  library(mixOmics)
  library(VennDiagram)
  library(grid)
  library(scales)
})

cat("\n=== FIGURE 2 - PHOSPHO LAYER ===\n")

# ----------------------------------------------------------------------------
# Restrict to the 3 regions that have phospho data (drop MockWholeBrain)
# ----------------------------------------------------------------------------
PHOS_TISSUES <- c("Cerebrum", "BrainStem", "Cerebellum")
phos_pep <- peptides %>%
  filter(tissue %in% PHOS_TISSUES, phospho)

# Phospho-protein × batch matrix (shared helper handles NA/zero filtering)
phos_mat <- build_protein_matrix(phos_pep)
phos_meta <- batch_meta %>% filter(tissue %in% PHOS_TISSUES)

cat(sprintf("Phospho matrix: %d batches x %d phosphoproteins\n",
            nrow(phos_mat), ncol(phos_mat) - 1))

# ----------------------------------------------------------------------------
# Panel A — PCA on phospho-protein intensities
# ----------------------------------------------------------------------------
phos_meta_grp <- phos_meta %>% mutate(group = tissue)
prep_A <- prep_plsda(phos_mat, meta = phos_meta_grp,
                     drop_batches = MV_BAD_BATCHES, top_var = 1500)
X_A <- prep_A$X
pca <- prcomp(X_A, scale. = TRUE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- as_tibble(pca$x[, 1:2], rownames = "batch") %>%
  left_join(phos_meta, by = "batch")

panel_A <- ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = sex)) +
  geom_point(size = 7, stroke = 1.4) +
  geom_text_repel(aes(label = sprintf("%s\n(%s)", tissue, tolower(as.character(sex)))),
                  size = 3.6, color = "grey20", fontface = "bold",
                  lineheight = 0.85,
                  box.padding = 0.7, point.padding = 0.5,
                  min.segment.length = 0.1, segment.color = "grey60",
                  max.overlaps = Inf) +
  scale_color_manual(values = TISSUE_COLORS, breaks = PHOS_TISSUES) +
  scale_shape_manual(values = c(Male = 16, Female = 17)) +
  labs(title = "PCA: phospho-proteome",
       subtitle = "Top 1,500 phosphoproteins by variance / Cerebellum-F: QC drop",
       x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2]),
       color = "Tissue", shape = "Sex") +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel B — PLS-DA with permutation validation (label: tissue)
# ----------------------------------------------------------------------------
set.seed(2026)
prep_B <- prep_plsda(phos_mat, meta = phos_meta_grp,
                     drop_batches = MV_BAD_BATCHES, top_var = 1500)
plsda_res <- plsda_with_perm(prep_B$X, prep_B$Y, n_perm = 999, ncomp = 2)

plsda_scores <- plsda_res$scores %>%
  left_join(phos_meta, by = "batch")

sig_label <- ifelse(plsda_res$p_val < 0.05, "validated", "n.s. (small n)")
panel_B <- ggplot(plsda_scores, aes(LV1, LV2, color = tissue, shape = sex)) +
  geom_point(size = 7, stroke = 1.4) +
  geom_text_repel(aes(label = sprintf("%s\n(%s)", tissue, tolower(as.character(sex)))),
                  size = 3.6, color = "grey20", fontface = "bold",
                  lineheight = 0.85,
                  box.padding = 0.7, point.padding = 0.5,
                  min.segment.length = 0.1, segment.color = "grey60",
                  max.overlaps = Inf) +
  scale_color_manual(values = TISSUE_COLORS, breaks = PHOS_TISSUES) +
  scale_shape_manual(values = c(Male = 16, Female = 17)) +
  labs(title = sprintf("PLS-DA: supervised on tissue (%s)", sig_label),
       subtitle = sprintf("999-perm permutation p = %.3f / n=5 batches",
                          plsda_res$p_val),
       x = sprintf("LV1 (%.1f%% var)", 100 * plsda_res$var_explained[1]),
       y = sprintf("LV2 (%.1f%% var)", 100 * plsda_res$var_explained[2]),
       color = "Tissue", shape = "Sex") +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel C — Top phospho-protein region markers (heatmap)
# ----------------------------------------------------------------------------
# One-vs-rest peptide-level Welch test on phospho subset, top 6 markers per region
markers <- bind_rows(lapply(PHOS_TISSUES, function(tt) {
  tissue_de_one_vs_rest(phos_pep, tt, channel = "heavy_intensity", min_n = 3)
}))

# Pick top 6 by adj.p (and require log2FC > 0 = up in target region)
top_markers <- markers %>%
  filter(!is.na(adj_p), log2FC > 0, !is.na(gene), gene != "") %>%
  group_by(tissue) %>%
  arrange(adj_p, desc(log2FC), .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup()

# Build z-scored intensity matrix on top markers
heat_data <- phos_pep %>%
  filter(protein_id %in% top_markers$protein_id, !is.na(heavy_intensity), heavy_intensity > 0) %>%
  mutate(log_int = log10(heavy_intensity)) %>%
  group_by(protein_id, gene, tissue, sex, batch) %>%
  summarise(value = mean(log_int, na.rm = TRUE), .groups = "drop") %>%
  group_by(protein_id) %>%
  mutate(z = scale(value)[, 1]) %>%
  ungroup() %>%
  inner_join(top_markers %>% dplyr::select(protein_id, marker_for = tissue),
             by = "protein_id") %>%
  mutate(
    sample = paste0(tissue, "\n", substr(as.character(sex), 1, 1)),
    sample = factor(sample,
                    levels = unique(paste0(rep(PHOS_TISSUES, each = 2),
                                           "\n",
                                           rep(c("M", "F"), times = 3))))
  )

heat_data$gene <- factor(
  heat_data$gene,
  levels = top_markers %>%
    arrange(factor(tissue, levels = rev(PHOS_TISSUES)), desc(log2FC)) %>%
    pull(gene) %>% unique()
)

panel_C <- ggplot(heat_data,
                  aes(x = sample, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2c7fb8", mid = "white", high = "#e31a1c",
                       midpoint = 0, na.value = "grey90",
                       name = "z-score\nlog10(heavy)") +
  facet_grid(rows = vars(marker_for), scales = "free_y", space = "free_y",
             switch = "y") +
  labs(title = "Top phospho-protein region markers",
       subtitle = "One-vs-rest Welch, top 6 per region by adj.p (log2FC > 0)",
       x = NULL, y = NULL) +
  theme_poster() +
  theme(
    axis.text.y       = element_text(face = "italic"),
    axis.text.x       = element_text(size = 10),
    strip.background  = element_rect(fill = "grey95", color = "grey70"),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 90, hjust = 0.5, face = "bold",
                                     size = 11, margin = margin(0, 4, 0, 4)),
    panel.spacing.y   = unit(4, "pt"),
    legend.position   = "right",
    legend.key.height = unit(1.0, "cm")
  )

# ----------------------------------------------------------------------------
# Panel D — Bulk vs phospho coverage Venn (per region)
# ----------------------------------------------------------------------------
bulk_prots_by_region <- peptides %>%
  filter(tissue %in% PHOS_TISSUES) %>%
  distinct(tissue, protein_id)
phos_prots_by_region <- phos_pep %>% distinct(tissue, protein_id)

cov_df <- map_dfr(PHOS_TISSUES, function(tt) {
  bulk <- bulk_prots_by_region %>% filter(tissue == tt) %>% pull(protein_id)
  phos <- phos_prots_by_region %>% filter(tissue == tt) %>% pull(protein_id)
  tibble(
    tissue = tt,
    bulk_only = length(setdiff(bulk, phos)),
    bulk_phos = length(intersect(bulk, phos)),
    phos_only = length(setdiff(phos, bulk))
  )
})
print(cov_df)

# "Unmod. only"  = proteins observed only via non-phospho peptides
# "Both layers"  = proteins observed in BOTH non-phospho and phospho peptides
# (Phospho-only is empirically 0 in this dataset and is dropped.)
cov_long <- cov_df %>%
  pivot_longer(c(bulk_only, bulk_phos),
               names_to = "category", values_to = "n") %>%
  mutate(
    category = recode(category,
                      bulk_only = "Unmodified only",
                      bulk_phos = "Both layers"),
    category = factor(category, levels = c("Unmodified only", "Both layers")),
    tissue = factor(tissue, levels = PHOS_TISSUES)
  )

cov_palette <- c("Unmodified only" = "#9ecae1",
                 "Both layers"     = "#3182bd")

panel_D <- ggplot(cov_long, aes(x = tissue, y = n, fill = category)) +
  geom_col(position = position_stack(reverse = TRUE),
           color = "white", linewidth = 0.4, width = 0.72) +
  geom_text(aes(label = format(n, big.mark = ",")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            color = "white", fontface = "bold", size = 3.6) +
  scale_fill_manual(values = cov_palette, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Phospho vs unmodified-peptide coverage per region",
       subtitle = "Unmodified-only vs both-layers proteins (phospho-only = 0 in this dataset)",
       x = NULL, y = "Proteins") +
  theme_poster() +
  theme(legend.position    = "top",
        legend.box.spacing = unit(2, "pt"),
        legend.spacing.x   = unit(14, "pt"),
        legend.key.width   = unit(14, "pt"),
        legend.text        = element_text(size = 11, margin = margin(l = 4, r = 12)))

# ----------------------------------------------------------------------------
# Combine (2 × 2)
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Phospho-proteome Layer",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("3 brain regions with phospho enrichment / Cerebrum, BrainStem, Cerebellum (Cerebellum-F: QC drop-out)",
             size = 13, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

top_row    <- plot_grid(panel_A, panel_B, ncol = 2,
                        labels = c("a", "b"), label_size = 18,
                        label_x = 0.01, label_y = 0.99)
bottom_row <- plot_grid(panel_C, panel_D, ncol = 2,
                        labels = c("c", "d"), label_size = 18,
                        label_x = 0.01, label_y = 0.99,
                        rel_widths = c(1.2, 1))

full <- plot_grid(header,
                  plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1.05)),
                  ncol = 1, rel_heights = c(0.07, 1))

save_poster(full, "FIG2_phospho", width = 16, height = 12)

# Persist tables for QA / supplementary
write_csv(top_markers,
          file.path(TAB_DIR, "FIG2_phospho_region_markers.csv"))
write_csv(cov_df,
          file.path(TAB_DIR, "FIG2_phospho_bulk_coverage.csv"))

cat(sprintf("PLS-DA permutation p = %.3f\n", plsda_res$p_val))
cat("\n=== FIGURE 2 written ===\n")
cat("  ../figures_poster/FIG2_phospho.pdf\n")
cat("  ../figures_poster/FIG2_phospho.png\n")
