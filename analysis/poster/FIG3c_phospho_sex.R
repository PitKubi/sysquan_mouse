#!/usr/bin/env Rscript
# FIGURE 3c — Sex differential on PHOSPHO-PROTEINS (CNPN 2026)
# Phosphopeptide-level intensities aggregated per protein × sex × peptide,
# then peptide-level Welch within protein. Heavy intensity is batch-median
# normalized (see 00_common.R). Restricted to Cerebrum + BrainStem for the
# volcanoes/concordance; Cerebellum phospho is included in the GO union.
#
# Panels:
#   a. Cerebrum  phospho-protein Male vs Female volcano
#   b. BrainStem phospho-protein Male vs Female volcano
#   c. Cross-region concordance scatter (log2FC Cer vs BS, phospho-proteins)
#   d. GO BP/MF/CC enrichment of phospho-sex genes (Cer + BS + Cbl union)
#
# A separate phosphoPEPTIDE-level figure (FIG3d) plots the same data at
# peptide resolution.
#
# Run from analysis/poster/:  Rscript FIG3c_phospho_sex.R

source("00_common.R")
suppressPackageStartupMessages({
  library(scales)
  library(ggrepel)
})

cat("\n=== FIGURE 3c - PHOSPHO SEX VOLCANOES + CONCORDANCE ===\n")

# Phospho subset
phos_pep <- peptides %>% filter(phospho)

# Per-region peptide-level sex DE on phospho. min_n=1 allows proteins with
# only one phosphopeptide per sex into the test (~5x more p-values computed)
# at the cost of unstable variance estimates — acceptable here as p-values
# are pseudo-replicate-driven and we use them as a ranking aid alongside
# log2FC rather than as a hard significance call.
sex_phos_cer <- sex_de_within_tissue(phos_pep, "Cerebrum",
                                     channel = "heavy_intensity",
                                     min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)
sex_phos_bs  <- sex_de_within_tissue(phos_pep, "BrainStem",
                                     channel = "heavy_intensity",
                                     min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)
sex_phos_cbl <- sex_de_within_tissue(phos_pep, "Cerebellum",
                                     channel = "heavy_intensity",
                                     min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)

tag_direction <- function(d, p_thresh = 0.05) {
  d %>%
    mutate(
      direction = case_when(
        is.na(p_value) | p_value >= p_thresh ~ "n.s.",
        log2FC > 0  ~ "Up in Male",
        log2FC < 0  ~ "Up in Female",
        TRUE        ~ "n.s."
      ),
      direction = factor(direction, levels = c("Up in Male", "Up in Female", "n.s."))
    )
}

cer_v <- tag_direction(sex_phos_cer)
bs_v  <- tag_direction(sex_phos_bs)
cbl_v <- tag_direction(sex_phos_cbl)

cat(sprintf("Phospho Cerebrum:   %d up-M / %d up-F / %d n.s.\n",
            sum(cer_v$direction == "Up in Male"),
            sum(cer_v$direction == "Up in Female"),
            sum(cer_v$direction == "n.s.")))
cat(sprintf("Phospho BrainStem:  %d up-M / %d up-F / %d n.s.\n",
            sum(bs_v$direction == "Up in Male"),
            sum(bs_v$direction == "Up in Female"),
            sum(bs_v$direction == "n.s.")))
cat(sprintf("Phospho Cerebellum: %d up-M / %d up-F / %d n.s.\n",
            sum(cbl_v$direction == "Up in Male"),
            sum(cbl_v$direction == "Up in Female"),
            sum(cbl_v$direction == "n.s.")))

# ----------------------------------------------------------------------------
# Volcano helper - protein-level (one point per phosphoprotein).
# Uses protein-level Welch test on batch-normalized log10(heavy).
# ----------------------------------------------------------------------------
volcano_panel <- function(d, region_name, n_label = 12) {
  d_plot <- d %>% filter(!is.na(log2FC), !is.na(p_value), p_value > 0)
  # Label top |log2FC| with reasonable signal (irrespective of p) — phospho
  # sex effects are weak after normalization, so labelling only p-sig
  # leaves the panel empty and uninformative.
  d_lab  <- d_plot %>%
    filter(!is.na(gene), gene != "", abs(log2FC) > 0.4) %>%
    arrange(desc(abs(log2FC))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = n_label)
  ggplot(d_plot, aes(x = log2FC, y = -log10(p_value), color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.55) +
    geom_point(data = d_lab, size = 2.3, alpha = 1) +
    geom_text_repel(data = d_lab, aes(label = gene),
                    size = 3.4, color = "grey15", fontface = "italic",
                    box.padding = 0.45, point.padding = 0.3,
                    min.segment.length = 0.1, segment.color = "grey55",
                    max.overlaps = Inf, seed = 1) +
    scale_color_manual(values = DIRECTION_COLORS, drop = FALSE) +
    labs(title = sprintf("%s phospho: Male vs Female", region_name),
         subtitle = sprintf("%d phosphoproteins / batch-normalized heavy / top |log2FC| labelled",
                            nrow(d_plot)),
         x = expression(log[2]~"(Male / Female)"),
         y = expression(-log[10]~"(p)"),
         color = NULL) +
    theme_poster() +
    theme(legend.position = "right")
}

panel_A <- volcano_panel(cer_v, "Cerebrum")
panel_B <- volcano_panel(bs_v,  "BrainStem")

# ----------------------------------------------------------------------------
# Panel C — Cross-region concordance scatter (Cer vs BS, phospho)
# ----------------------------------------------------------------------------
concord <- inner_join(
  sex_phos_cer %>% dplyr::select(protein_id, gene, log2FC_Cer = log2FC,
                                 p_Cer = p_value),
  sex_phos_bs  %>% dplyr::select(protein_id, log2FC_BS = log2FC,
                                 p_BS = p_value),
  by = "protein_id"
) %>%
  filter(!is.na(log2FC_Cer), !is.na(log2FC_BS), !is.na(gene), gene != "")

concord <- concord %>%
  mutate(
    sig_in_either = pmin(p_Cer, p_BS, na.rm = TRUE) < 0.05,
    concordant    = sig_in_either & (sign(log2FC_Cer) == sign(log2FC_BS)),
    class = case_when(
      concordant & log2FC_Cer > 0 ~ "Shared - Up in Male",
      concordant & log2FC_Cer < 0 ~ "Shared - Up in Female",
      sig_in_either               ~ "Region-specific",
      TRUE                        ~ "n.s."
    ),
    class = factor(class, levels = c("Shared - Up in Male",
                                     "Shared - Up in Female",
                                     "Region-specific", "n.s."))
  )

# Label top 10 concordant + 5 strongest off-diagonal
lab_concord <- concord %>% filter(concordant) %>%
  arrange(desc(abs((log2FC_Cer + log2FC_BS) / 2))) %>% slice_head(n = 10)
lab_offdiag <- concord %>% filter(class == "Region-specific") %>%
  arrange(desc(abs(log2FC_Cer - log2FC_BS))) %>% slice_head(n = 5)
lab_set <- bind_rows(lab_concord, lab_offdiag) %>% distinct(protein_id, .keep_all = TRUE)

class_pal <- c(
  "Shared - Up in Male"   = SEX_COLORS[["Male"]],
  "Shared - Up in Female" = SEX_COLORS[["Female"]],
  "Region-specific"       = "#fdae61",
  "n.s."                  = "grey80"
)

axis_lim <- max(abs(c(concord$log2FC_Cer, concord$log2FC_BS)), na.rm = TRUE)

panel_C <- ggplot(concord, aes(x = log2FC_Cer, y = log2FC_BS, color = class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.4) +
  geom_point(data = filter(concord, class == "n.s."), size = 1.0, alpha = 0.4) +
  geom_point(data = filter(concord, class != "n.s."), size = 2.2, alpha = 0.85) +
  geom_text_repel(data = lab_set, aes(label = gene),
                  size = 3.3, color = "grey15", fontface = "italic",
                  box.padding = 0.45, point.padding = 0.3,
                  min.segment.length = 0.1, segment.color = "grey55",
                  max.overlaps = Inf, seed = 2) +
  scale_color_manual(values = class_pal, drop = FALSE) +
  scale_x_continuous(limits = c(-axis_lim, axis_lim)) +
  scale_y_continuous(limits = c(-axis_lim, axis_lim)) +
  coord_fixed() +
  labs(title = "Cross-region concordance of phospho sex effects",
       subtitle = "Diagonal = identical effect in both regions; sig. = peptide-level p<0.05 in at least one",
       x = expression(log[2]~"(M/F) - Cerebrum (phospho)"),
       y = expression(log[2]~"(M/F) - BrainStem (phospho)"),
       color = NULL) +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel D — GO enrichment (BP + MF + CC) of phospho sex genes
# Combine two selection criteria per region (union):
#   - phosphoprotein-level Welch p-value < 0.30 (after relaxing min_n=1), OR
#   - top 400 by |log2FC| per region (always-available ranking signal).
# This intentionally produces a larger input list than a hard FDR/p threshold
# would; clusterProfiler then applies its own BH-FDR over GO terms.
# ----------------------------------------------------------------------------
loose_sig_set <- function(d, p_thresh = 0.30, top_n = 400) {
  by_p  <- d %>% filter(!is.na(p_value), p_value < p_thresh,
                        !is.na(log2FC), is.finite(log2FC))
  by_fc <- d %>% filter(!is.na(log2FC), is.finite(log2FC)) %>%
    arrange(desc(abs(log2FC))) %>% slice_head(n = top_n)
  bind_rows(by_p, by_fc) %>% distinct(protein_id, .keep_all = TRUE)
}

sig_union <- bind_rows(
  loose_sig_set(sex_phos_cer) %>% mutate(region = "Cerebrum"),
  loose_sig_set(sex_phos_bs)  %>% mutate(region = "BrainStem"),
  loose_sig_set(sex_phos_cbl) %>% mutate(region = "Cerebellum")
)

sex_genes_by_dir <- sig_union %>%
  filter(!is.na(gene), gene != "") %>%
  group_by(gene) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(direction = ifelse(mean_log2FC > 0, "Up in Male", "Up in Female"))

cat(sprintf("Phospho sex genes (loose): %d up-M / %d up-F (Cer+BS+Cbl union)\n",
            sum(sex_genes_by_dir$direction == "Up in Male"),
            sum(sex_genes_by_dir$direction == "Up in Female")))

# Phospho-aware universe
phos_universe <- phos_pep %>%
  filter(!is.na(gene), gene != "") %>% pull(gene) %>% unique()
cat(sprintf("Phospho universe: %d genes\n", length(phos_universe)))

go_sex <- bind_rows(lapply(c("Up in Male", "Up in Female"), function(dir) {
  syms <- sex_genes_by_dir %>% filter(direction == dir) %>% pull(gene) %>% unique()
  enrich_gene_set_all(syms, label = dir, fdr_cap = 0.25, top_n = 5,
                      universe = phos_universe)
}))

panel_D <- plot_go_3panels(
  go_sex,
  title    = "GO enrichment - phospho sex-effect genes",
  subtitle = sprintf("BP / MF / CC top 5 each; %d up-M, %d up-F genes (Cer+BS+Cbl union, p<0.30 OR top 400 |log2FC|)",
                     sum(sex_genes_by_dir$direction == "Up in Male"),
                     sum(sex_genes_by_dir$direction == "Up in Female"))
)

# ----------------------------------------------------------------------------
# Combine: 2 × 2 grid
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Sex Effect on the Phospho-proteome",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("Cerebrum + BrainStem volcanoes and concordance / GO enrichment includes Cerebellum",
             size = 13, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

row1 <- plot_grid(panel_A, panel_B, ncol = 2,
                  labels = c("a", "b"), label_size = 18,
                  label_x = 0.01, label_y = 0.99)
row2 <- plot_grid(panel_C, panel_D, ncol = 2,
                  labels = c("c", "d"), label_size = 18,
                  label_x = 0.01, label_y = 0.99,
                  rel_widths = c(1, 1.05))

full <- plot_grid(header,
                  plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.05)),
                  ncol = 1, rel_heights = c(0.05, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_poster(full, "FIG3c_phospho_sex", width = 17, height = 17)

write_csv(concord %>% arrange(desc(concordant), p_Cer + p_BS),
          file.path(TAB_DIR, "FIG3c_phospho_sex_concordance.csv"))
write_csv(go_sex, file.path(TAB_DIR, "FIG3c_GO_phospho_sex.csv"))

cat(sprintf("Concordant phospho proteins: %d\n", sum(concord$concordant, na.rm = TRUE)))
cat("\n=== FIGURE 3c written ===\n")
cat(sprintf("  %s/FIG3c_phospho_sex.pdf\n", FIG_DIR))
cat(sprintf("  %s/FIG3c_phospho_sex.png\n", FIG_DIR))
