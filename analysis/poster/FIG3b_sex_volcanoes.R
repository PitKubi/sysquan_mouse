#!/usr/bin/env Rscript
# FIGURE 3b — Sex differential deep-dive on the PROTEOME (CNPN 2026)
# Restricted to Cerebrum + BrainStem (only paired regions for sex).
# Panels:
#   a. Cerebrum  Male vs Female volcano (peptide-level Welch within protein)
#   b. BrainStem Male vs Female volcano
#   c. Cross-tissue concordance scatter (log2FC Cer vs BS)
#   d. GO BP enrichment of significant sex genes, split by direction
#
# Run from analysis/poster/:  Rscript FIG3b_sex_volcanoes.R

source("00_common.R")
suppressPackageStartupMessages({
  library(scales)
  library(ggrepel)
})

cat("\n=== FIGURE 3b - SEX VOLCANOES + CONCORDANCE ===\n")

# ----------------------------------------------------------------------------
# Per-region peptide-level sex DE
# ----------------------------------------------------------------------------
sex_de_cer <- sex_de_within_tissue(peptides, "Cerebrum",  channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)
sex_de_bs  <- sex_de_within_tissue(peptides, "BrainStem", channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)

# Tag direction at p<0.05 (peptide-level, no FDR — honest given small n)
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

cer_v <- tag_direction(sex_de_cer)
bs_v  <- tag_direction(sex_de_bs)

cat(sprintf("Cerebrum:  %d up-M / %d up-F / %d n.s.  (peptide-Welch p<0.05)\n",
            sum(cer_v$direction == "Up in Male"),
            sum(cer_v$direction == "Up in Female"),
            sum(cer_v$direction == "n.s.")))
cat(sprintf("BrainStem: %d up-M / %d up-F / %d n.s.\n",
            sum(bs_v$direction == "Up in Male"),
            sum(bs_v$direction == "Up in Female"),
            sum(bs_v$direction == "n.s.")))

# ----------------------------------------------------------------------------
# Volcano helper
# ----------------------------------------------------------------------------
volcano_panel <- function(d, region_name, n_label = 10) {
  d_plot <- d %>% filter(!is.na(log2FC), !is.na(p_value), p_value > 0)
  d_lab  <- d_plot %>%
    filter(direction != "n.s.", !is.na(gene), gene != "") %>%
    arrange(p_value, desc(abs(log2FC))) %>%
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
    labs(title = sprintf("%s: Male vs Female", region_name),
         subtitle = sprintf("%d proteins tested / peptide-level Welch within protein",
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
# Panel C — Cross-tissue concordance scatter
# ----------------------------------------------------------------------------
concord <- inner_join(
  sex_de_cer %>% dplyr::select(protein_id, gene, log2FC_Cer = log2FC,
                               p_Cer = p_value, adj_p_Cer = adj_p),
  sex_de_bs  %>% dplyr::select(protein_id, log2FC_BS = log2FC,
                               p_BS = p_value, adj_p_BS = adj_p),
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

# Label top 12 concordant + 5 strongest off-diagonal
lab_concord <- concord %>%
  filter(concordant) %>%
  arrange(desc(abs((log2FC_Cer + log2FC_BS) / 2))) %>%
  slice_head(n = 12)
lab_offdiag <- concord %>%
  filter(class == "Region-specific") %>%
  arrange(desc(abs(log2FC_Cer - log2FC_BS))) %>%
  slice_head(n = 5)
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
  geom_point(data = filter(concord, class == "n.s."),
             size = 1.0, alpha = 0.4) +
  geom_point(data = filter(concord, class != "n.s."),
             size = 2.2, alpha = 0.85) +
  geom_text_repel(data = lab_set, aes(label = gene),
                  size = 3.3, color = "grey15", fontface = "italic",
                  box.padding = 0.45, point.padding = 0.3,
                  min.segment.length = 0.1, segment.color = "grey55",
                  max.overlaps = Inf, seed = 2) +
  scale_color_manual(values = class_pal, drop = FALSE) +
  scale_x_continuous(limits = c(-axis_lim, axis_lim)) +
  scale_y_continuous(limits = c(-axis_lim, axis_lim)) +
  coord_fixed() +
  labs(title = "Cross-region concordance of sex effects",
       subtitle = "Diagonal = identical effect in both regions; sig. = peptide-level p<0.05 in at least one region",
       x = expression(log[2]~"(M/F) - Cerebrum"),
       y = expression(log[2]~"(M/F) - BrainStem"),
       color = NULL) +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel D — GO enrichment (BP + MF + CC) of sex-effect genes (Cer + BS union)
# Loose cutoff: peptide-level p<0.10 OR top 200 by |log2FC| per region.
# ----------------------------------------------------------------------------
loose_sig_set <- function(d, p_thresh = 0.10, top_n = 200) {
  by_p <- d %>% filter(!is.na(p_value), p_value < p_thresh, !is.na(log2FC))
  by_fc <- d %>% filter(!is.na(log2FC)) %>%
    arrange(desc(abs(log2FC))) %>% slice_head(n = top_n)
  bind_rows(by_p, by_fc) %>% distinct(protein_id, .keep_all = TRUE)
}

sig_union <- bind_rows(
  loose_sig_set(sex_de_cer) %>% mutate(region = "Cerebrum"),
  loose_sig_set(sex_de_bs)  %>% mutate(region = "BrainStem")
)

sex_genes_by_dir <- sig_union %>%
  filter(!is.na(gene), gene != "") %>%
  group_by(gene) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(direction = ifelse(mean_log2FC > 0, "Up in Male", "Up in Female"))

cat(sprintf("Sex genes (proteome, loose): %d up-M / %d up-F (Cer+BS union)\n",
            sum(sex_genes_by_dir$direction == "Up in Male"),
            sum(sex_genes_by_dir$direction == "Up in Female")))

go_sex <- bind_rows(lapply(c("Up in Male", "Up in Female"), function(dir) {
  syms <- sex_genes_by_dir %>% filter(direction == dir) %>% pull(gene) %>% unique()
  enrich_gene_set_all(syms, label = dir, fdr_cap = 0.25, top_n = 5)
}))

panel_D <- plot_go_3panels(
  go_sex,
  title    = "GO enrichment - sex-effect genes (proteome)",
  subtitle = sprintf("BP / MF / CC top 5 each; %d up-M, %d up-F genes (Cer+BS union, p<0.10 OR top 200 |log2FC|)",
                     sum(sex_genes_by_dir$direction == "Up in Male"),
                     sum(sex_genes_by_dir$direction == "Up in Female"))
)

# ----------------------------------------------------------------------------
# Combine: 2 × 2 grid (a | b) / (c | d)
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Sex Effect on the Proteome (Cerebrum + BrainStem)",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("Per-region volcanoes, cross-region concordance, and GO enrichment of M-vs-F log2FC",
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

save_poster(full, "FIG3b_sex_volcanoes", width = 17, height = 17)
write_csv(go_sex, file.path(TAB_DIR, "FIG3b_GO_sex_proteome.csv"))

write_csv(concord %>%
            arrange(desc(concordant), p_Cer + p_BS),
          file.path(TAB_DIR, "FIG3b_sex_concordance_full.csv"))

cat(sprintf("Concordant (sig+sign-match) proteins: %d\n",
            sum(concord$concordant, na.rm = TRUE)))
cat(sprintf("Region-specific sig: %d\n",
            sum(concord$class == "Region-specific", na.rm = TRUE)))
cat("\n=== FIGURE 3b written ===\n")
cat(sprintf("  %s/FIG3b_sex_volcanoes.pdf\n", FIG_DIR))
cat(sprintf("  %s/FIG3b_sex_volcanoes.png\n", FIG_DIR))
