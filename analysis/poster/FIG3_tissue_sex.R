#!/usr/bin/env Rscript
# FIGURE 3 — Tissue and sex differential drivers (CNPN 2026)
# Panels:
#   a. PCA of full proteome (6 healthy batches), tissue × sex
#   b. PLS-DA on full proteome, tissue label, with permutation p-value
#   c. UpSet of tissue marker proteins (top 200/region one-vs-rest)
#   d. Tissue driver heatmap (top 5 markers per region, z-scored)
#   e. Sex driver heatmap (concordant log2FC across Cerebrum + BrainStem)
#
# Run from analysis/poster/:  Rscript FIG3_tissue_sex.R

source("00_common.R")
suppressPackageStartupMessages({
  library(mixOmics)
  library(ComplexUpset)
  library(scales)
})

cat("\n=== FIGURE 3 - TISSUE x SEX DIFFERENTIAL ===\n")

# ----------------------------------------------------------------------------
# Use full proteome (all 4 tissues, drop QC-fail batch for multivariate)
# ----------------------------------------------------------------------------
full_meta <- batch_meta %>% mutate(group = tissue)

# ----------------------------------------------------------------------------
# Panel A — PCA, full proteome
# ----------------------------------------------------------------------------
prep_full <- prep_plsda(batch_matrix_full, meta = full_meta,
                        drop_batches = MV_BAD_BATCHES, top_var = 2000)
pca <- prcomp(prep_full$X, scale. = TRUE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- as_tibble(pca$x[, 1:2], rownames = "batch") %>%
  left_join(full_meta, by = "batch")

panel_A <- ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = sex)) +
  geom_point(size = 7, stroke = 1.4) +
  geom_text_repel(aes(label = sprintf("%s\n(%s)", tissue, tolower(as.character(sex)))),
                  size = 3.6, color = "grey20", fontface = "bold",
                  lineheight = 0.85,
                  box.padding = 0.7, point.padding = 0.5,
                  min.segment.length = 0.1, segment.color = "grey60",
                  max.overlaps = Inf) +
  scale_color_manual(values = TISSUE_COLORS) +
  scale_shape_manual(values = c(Male = 16, Female = 17)) +
  labs(title = "PCA: full proteome",
       subtitle = "Top 2,000 proteins by variance / Cerebellum-F: QC drop",
       x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2]),
       color = "Tissue", shape = "Sex") +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel B — PLS-DA on tissue, permutation-validated
# ----------------------------------------------------------------------------
set.seed(2026)
plsda_res <- plsda_with_perm(prep_full$X, prep_full$Y, n_perm = 999, ncomp = 2)
plsda_scores <- plsda_res$scores %>% left_join(full_meta, by = "batch")

sig_label <- ifelse(plsda_res$p_val < 0.05, "validated", "n.s. (small n)")
panel_B <- ggplot(plsda_scores, aes(LV1, LV2, color = tissue, shape = sex)) +
  geom_point(size = 7, stroke = 1.4) +
  geom_text_repel(aes(label = sprintf("%s\n(%s)", tissue, tolower(as.character(sex)))),
                  size = 3.6, color = "grey20", fontface = "bold",
                  lineheight = 0.85,
                  box.padding = 0.7, point.padding = 0.5,
                  min.segment.length = 0.1, segment.color = "grey60",
                  max.overlaps = Inf) +
  scale_color_manual(values = TISSUE_COLORS) +
  scale_shape_manual(values = c(Male = 16, Female = 17)) +
  labs(title = sprintf("PLS-DA: supervised on tissue (%s)", sig_label),
       subtitle = sprintf("999-perm p = %.3f / n=6 batches",
                          plsda_res$p_val),
       x = sprintf("LV1 (%.1f%% var)", 100 * plsda_res$var_explained[1]),
       y = sprintf("LV2 (%.1f%% var)", 100 * plsda_res$var_explained[2]),
       color = "Tissue", shape = "Sex") +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel C — UpSet of tissue marker proteins (top 200/region one-vs-rest)
# ----------------------------------------------------------------------------
tissue_de <- bind_rows(lapply(levels(peptides$tissue), function(tt) {
  tissue_de_one_vs_rest(peptides, tt, channel = "heavy_intensity", min_n = 3)
}))

top_tissue_markers <- tissue_de %>%
  filter(!is.na(adj_p), log2FC > 0, !is.na(gene), gene != "") %>%
  group_by(tissue) %>%
  arrange(adj_p, desc(log2FC), .by_group = TRUE) %>%
  slice_head(n = 200) %>%
  ungroup()

# Wide presence matrix
marker_wide <- top_tissue_markers %>%
  distinct(tissue, protein_id) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = tissue, values_from = present, values_fill = FALSE) %>%
  mutate(across(-protein_id, as.logical))

ups_tissues <- intersect(names(TISSUE_COLORS), colnames(marker_wide))

panel_C <- ComplexUpset::upset(
  marker_wide,
  intersect = ups_tissues,
  name = "",
  width_ratio = 0.20, height_ratio = 0.55,
  min_size = 1, n_intersections = 12,
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Markers" = intersection_size(
      text = list(size = 3.2, vjust = -0.2),
      mapping = aes(fill = "bar_color")
    ) +
      scale_fill_manual(values = c(bar_color = "#33415C"), guide = "none") +
      labs(y = "Markers in intersection") +
      theme_poster(base_size = 12) +
      theme(panel.grid.major.x = element_blank())
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(aes(fill = group), color = NA, width = 0.7)
  ) +
    scale_fill_manual(values = TISSUE_COLORS, guide = "none") +
    labs(x = "Per-tissue markers (top 200)") +
    theme_poster(base_size = 11),
  matrix = intersection_matrix(
    geom = geom_point(size = 2.3),
    segment = geom_segment(linewidth = 0.7)
  ),
  themes = upset_default_themes(text = element_text(size = 12))
) +
  plot_annotation(title = "Region-marker overlap (top 200/region, one-vs-rest)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

# ----------------------------------------------------------------------------
# Panel D — Tissue driver heatmap (top 5 markers per region)
# ----------------------------------------------------------------------------
top5 <- top_tissue_markers %>%
  group_by(tissue) %>% slice_head(n = 4) %>% ungroup() %>%
  distinct(protein_id, .keep_all = TRUE)  # if a protein is top in 2 regions, keep first


tissue_heat <- peptides %>%
  filter(protein_id %in% top5$protein_id, !is.na(heavy_intensity), heavy_intensity > 0,
         !batch %in% MV_BAD_BATCHES, !is.na(sex)) %>%
  mutate(log_int = log10(heavy_intensity)) %>%
  group_by(protein_id, gene, batch, tissue, sex) %>%
  summarise(value = mean(log_int, na.rm = TRUE), .groups = "drop") %>%
  group_by(protein_id) %>%
  mutate(z = scale(value)[, 1]) %>%
  ungroup() %>%
  inner_join(top5 %>% dplyr::select(protein_id, marker_for = tissue),
             by = "protein_id") %>%
  mutate(
    sample = paste0(tissue, "\n", substr(as.character(sex), 1, 1)),
    sample = factor(sample,
                    levels = c("Cerebrum\nM","Cerebrum\nF",
                               "BrainStem\nM","BrainStem\nF",
                               "Cerebellum\nM",
                               "MockWholeBrain\nM"))
  )

# Make gene labels unique (some gene symbols repeat across protein_ids)
tissue_heat <- tissue_heat %>%
  mutate(gene_label = make.unique(as.character(gene)))
ord <- top5 %>%
  arrange(factor(tissue, levels = rev(levels(peptides$tissue))), desc(log2FC)) %>%
  pull(protein_id)
tissue_heat$gene_label <- factor(tissue_heat$gene_label,
                                 levels = unique(tissue_heat$gene_label[order(match(tissue_heat$protein_id, ord))]))

# Use shorter tissue tags for the strip so they fit at angle=0
tissue_heat <- tissue_heat %>%
  mutate(marker_for = factor(recode(as.character(marker_for),
                                    Cerebrum = "Cerebrum\nmarker",
                                    BrainStem = "BrainStem\nmarker",
                                    Cerebellum = "Cerebellum\nmarker",
                                    MockWholeBrain = "Mock\nmarker"),
                             levels = c("Cerebrum\nmarker", "BrainStem\nmarker",
                                        "Cerebellum\nmarker", "Mock\nmarker")))

panel_D <- ggplot(tissue_heat, aes(x = sample, y = gene_label, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2c7fb8", mid = "white", high = "#e31a1c",
                       midpoint = 0, na.value = "grey90",
                       name = "z-score\nlog10(heavy)") +
  facet_grid(rows = vars(marker_for), scales = "free_y", space = "free_y",
             switch = "y") +
  labs(title = "Top 4 region markers per tissue",
       subtitle = "Heavy-channel intensity (z-scored across batches)",
       x = NULL, y = NULL) +
  theme_poster() +
  theme(
    axis.text.y       = element_text(face = "italic", size = 10),
    axis.text.x       = element_text(size = 10),
    strip.background  = element_rect(fill = "grey95", color = "grey70"),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold",
                                     size = 11, lineheight = 0.85,
                                     margin = margin(0, 6, 0, 4)),
    panel.spacing.y   = unit(4, "pt"),
    legend.position   = "right",
    legend.key.height = unit(1.0, "cm")
  )

# ----------------------------------------------------------------------------
# Panel E — Sex driver heatmap (concordant log2FC across Cerebrum + BrainStem)
# ----------------------------------------------------------------------------
sex_de_cer <- sex_de_within_tissue(peptides, "Cerebrum",  channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)
sex_de_bs  <- sex_de_within_tissue(peptides, "BrainStem", channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)

sex_concord <- inner_join(
  sex_de_cer %>% dplyr::select(protein_id, gene, log2FC_Cer = log2FC,
                               adj_p_Cer = adj_p),
  sex_de_bs  %>% dplyr::select(protein_id, log2FC_BS = log2FC, adj_p_BS = adj_p),
  by = "protein_id"
) %>%
  filter(!is.na(log2FC_Cer), !is.na(log2FC_BS),
         !is.na(gene), gene != "",
         sign(log2FC_Cer) == sign(log2FC_BS),
         pmin(adj_p_Cer, adj_p_BS, na.rm = TRUE) < 0.10) %>%
  mutate(
    mean_log2FC = (log2FC_Cer + log2FC_BS) / 2,
    direction   = ifelse(mean_log2FC > 0, "Up in Male", "Up in Female")
  )

# Top 8 per direction by |mean log2FC|
top_sex_concord <- sex_concord %>%
  group_by(direction) %>%
  arrange(desc(abs(mean_log2FC)), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

# Reshape for tile plot: rows = gene, cols = (Cerebrum, BrainStem), value = log2FC
sex_heat <- top_sex_concord %>%
  dplyr::select(gene, direction, log2FC_Cer, log2FC_BS) %>%
  pivot_longer(c(log2FC_Cer, log2FC_BS),
               names_to = "tissue", values_to = "log2FC") %>%
  mutate(tissue = recode(tissue, log2FC_Cer = "Cerebrum", log2FC_BS = "BrainStem"),
         tissue = factor(tissue, levels = c("Cerebrum", "BrainStem")),
         direction = recode(direction,
                            "Up in Male"   = "Up in\nMale",
                            "Up in Female" = "Up in\nFemale"),
         direction = factor(direction, levels = c("Up in\nMale", "Up in\nFemale")))

# Uniqueify any duplicate gene labels and lock factor order
sex_heat <- sex_heat %>%
  group_by(gene) %>%
  mutate(gene_label = if (n_distinct(direction) > 1)
                        paste0(gene, "_", as.integer(direction)) else as.character(gene)) %>%
  ungroup()
gene_order <- top_sex_concord %>%
  arrange(direction, desc(abs(mean_log2FC))) %>% pull(gene) %>% as.character()
sex_heat$gene_label <- factor(sex_heat$gene_label,
                              levels = rev(unique(sex_heat$gene_label[
                                order(match(sex_heat$gene, gene_order))])))

max_abs <- max(abs(sex_heat$log2FC), na.rm = TRUE)

panel_E <- ggplot(sex_heat, aes(x = tissue, y = gene_label, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.1f", log2FC)), size = 3.0,
            color = "grey15", fontface = "bold") +
  scale_fill_gradient2(low = SEX_COLORS["Female"], mid = "white",
                       high = SEX_COLORS["Male"], midpoint = 0,
                       limits = c(-max_abs, max_abs),
                       name = "log2(M/F)") +
  facet_grid(rows = vars(direction), scales = "free_y", space = "free_y",
             switch = "y") +
  labs(title = "Sex-concordant proteins (Cerebrum + BrainStem)",
       subtitle = "Top 8 per direction by mean |log2FC|; FDR < 0.10 in at least one region",
       x = NULL, y = NULL) +
  theme_poster() +
  theme(
    axis.text.y       = element_text(face = "italic", size = 10),
    strip.background  = element_rect(fill = "grey95", color = "grey70"),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold",
                                     size = 11, lineheight = 0.85,
                                     margin = margin(0, 6, 0, 4)),
    panel.spacing.y   = unit(4, "pt"),
    legend.position   = "right",
    legend.key.height = unit(1.0, "cm")
  )

# ----------------------------------------------------------------------------
# Combine: 2 rows top (A, B), 1 row middle (C UpSet half + E sex heat half)
# Tissue marker heatmap (panel D in earlier draft) was redundant with c/UpSet
# and GO-tissue panel in FIG4. Dropped to keep this figure clean.
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Tissue and Sex Differential",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("Drivers of regional and sex-dimorphic protein expression",
             size = 13, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

panel_C_grob <- cowplot::ggdraw(panel_C)

row1 <- plot_grid(panel_A, panel_B, ncol = 2,
                  labels = c("a", "b"), label_size = 18,
                  label_x = 0.01, label_y = 0.99)
row2 <- plot_grid(panel_C_grob, panel_E, ncol = 2,
                  labels = c("c", "d"), label_size = 18,
                  label_x = 0.01, label_y = 0.99,
                  rel_widths = c(1.15, 1))

full <- plot_grid(header,
                  plot_grid(row1, row2, ncol = 1,
                            rel_heights = c(1, 1.15)),
                  ncol = 1, rel_heights = c(0.05, 1))

save_poster(full, "FIG3_tissue_sex", width = 17, height = 13)

# Persist tables
write_csv(top_tissue_markers,
          file.path(TAB_DIR, "FIG3_tissue_markers_top200.csv"))
write_csv(top_sex_concord,
          file.path(TAB_DIR, "FIG3_sex_concordant_topN.csv"))

cat(sprintf("Tissue PLS-DA permutation p = %.3f\n", plsda_res$p_val))
cat(sprintf("Sex-concordant proteins (FDR<0.10): %d total / %d top selected\n",
            nrow(sex_concord), nrow(top_sex_concord)))
cat("\n=== FIGURE 3 written ===\n")
cat("  ../figures_poster/FIG3_tissue_sex.pdf\n")
cat("  ../figures_poster/FIG3_tissue_sex.png\n")
