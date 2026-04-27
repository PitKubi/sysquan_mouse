#!/usr/bin/env Rscript
# FIGURE 1 — Mouse Brain SysQuan: resource overview (CNPN 2026)
# Panels:
#   A. Pyramid funnel of peptide-level filtering (Total → Proteotypic →
#      Quantifiable doublet → Phospho-detectable)
#   B. Unique proteins per tissue × sex (Cerebellum_F flagged as QC drop-out)
#   C. UpSet plot of proteins detected across the 4 tissues
#   D. Peptides per protein distribution by tissue (coverage depth)
#
# Run from analysis/poster/:  Rscript FIG1_resource.R

source("00_common.R")
suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(scales)
})

cat("\n=== FIGURE 1 — RESOURCE OVERVIEW ===\n")

# ----------------------------------------------------------------------------
# Panel A — Protein-level pyramid (Total → Proteotypic → Quantifiable → Phospho)
# Each tier is a strict subset of the previous (computed at the protein level).
# ----------------------------------------------------------------------------
prot_total       <- peptides %>% pull(protein_id) %>% n_distinct()
prot_proteotypic <- peptides %>% filter(proteotypic) %>% pull(protein_id) %>% n_distinct()
prot_quant       <- peptides %>% filter(proteotypic, doublet) %>% pull(protein_id) %>% n_distinct()
phos_quant_set   <- intersect(
  peptides %>% filter(phospho) %>% pull(protein_id),
  peptides %>% filter(proteotypic, doublet) %>% pull(protein_id)
)
prot_phos_quant  <- length(unique(phos_quant_set))

# Peptide-level totals reported in the subtitle (context for clients reading)
pep_total <- nrow(peptides)
pep_phos  <- sum(peptides$phospho, na.rm = TRUE)

pyramid_df <- tibble(
  level = c("Total proteins detected",
            "Proteotypic",
            "Quantifiable",
            "Phospho-resolved"),
  sub   = c("all peptide observations",
            ">=1 proteotypic peptide",
            "proteotypic + doublet",
            "quantifiable + phosphosite"),
  count = c(prot_total, prot_proteotypic, prot_quant, prot_phos_quant)
) %>%
  mutate(
    y = rev(seq_along(level)),  # 4 = top, 1 = bottom
    pct_of_total = 100 * count / prot_total,
    half = (count / max(count)) / 2,
    xmin = 0.5 - half,
    xmax = 0.5 + half,
    is_wide = (xmax - xmin) > 0.55,
    level = factor(level, levels = rev(level))
  )

# Build connecting trapezoids for a funnel look
trap_rows <- pyramid_df %>%
  arrange(desc(y)) %>%
  mutate(next_xmin = lead(xmin), next_xmax = lead(xmax),
         next_y    = lead(y)) %>%
  filter(!is.na(next_y))

trap_poly <- trap_rows %>%
  rowwise() %>%
  do({
    r <- .
    tibble(
      group = paste0("trap_", r$y),
      poly_x = c(r$xmin, r$xmax, r$next_xmax, r$next_xmin),
      poly_y = c(r$y - 0.40, r$y - 0.40, r$next_y + 0.40, r$next_y + 0.40)
    )
  }) %>%
  ungroup()

pyramid_pal <- c("#264653", "#2A9D8F", "#E9C46A", "#E76F51")

panel_A <- ggplot() +
  # trapezoid connectors
  geom_polygon(data = trap_poly,
               aes(x = poly_x, y = poly_y, group = group),
               fill = "grey85", color = NA) +
  # rectangles for each tier
  geom_rect(data = pyramid_df,
            aes(xmin = xmin, xmax = xmax,
                ymin = y - 0.40, ymax = y + 0.40, fill = level),
            color = "white", linewidth = 0.6) +
  # tier name (always LEFT of plot, hjust=1)
  geom_text(data = pyramid_df,
            aes(x = -0.04, y = y + 0.08, label = level),
            hjust = 1, fontface = "bold", size = 4.5, color = "grey15") +
  geom_text(data = pyramid_df,
            aes(x = -0.04, y = y - 0.13, label = sub),
            hjust = 1, fontface = "italic", size = 3.3, color = "grey45") +
  # count + percentage (always RIGHT of plot, hjust=0)
  geom_text(data = pyramid_df,
            aes(x = 1.04, y = y + 0.04,
                label = format(count, big.mark = ",")),
            hjust = 0, fontface = "bold", size = 4.7, color = "grey15") +
  geom_text(data = pyramid_df,
            aes(x = 1.04, y = y - 0.18,
                label = sprintf("(%.1f%%)", pct_of_total)),
            hjust = 0, size = 3.5, color = "grey45") +
  scale_fill_manual(values = setNames(pyramid_pal, levels(pyramid_df$level))) +
  scale_y_continuous(limits = c(0.40, 4.60), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.95, 1.50), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Protein-level resource composition",
       subtitle = sprintf("%s peptide entries; %s phosphopeptides; 9,256 proteins shared across all 4 regions",
                          format(pep_total, big.mark = ","),
                          format(pep_phos,  big.mark = ","))) +
  theme_void(base_size = 13) +
  theme(legend.position = "none",
        plot.title       = element_text(size = 15, face = "bold", hjust = 0,
                                        margin = margin(0, 0, 2, 22)),
        plot.subtitle    = element_text(size = 10.5, color = "grey30", hjust = 0,
                                        margin = margin(0, 0, 6, 22)),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin      = margin(10, 10, 10, 10))

# ----------------------------------------------------------------------------
# Panel B — Unique proteins per tissue × sex (with QC drop-out flag)
# ----------------------------------------------------------------------------
B_data <- peptides %>%
  filter(!is.na(sex)) %>%
  group_by(tissue, sex, batch) %>%
  summarise(n_proteins = n_distinct(protein_id), .groups = "drop") %>%
  mutate(qc_fail = batch %in% MV_BAD_BATCHES)

# Annotation segment: "no fractionation" arrow pointing to Cerebellum-F bar
qc_row <- B_data %>% filter(qc_fail)

panel_B <- ggplot(B_data, aes(x = tissue, y = n_proteins, fill = sex)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = format(n_proteins, big.mark = ",")),
            position = position_dodge(width = 0.85), vjust = -0.4,
            size = 3.4, fontface = "bold") +
  # "no fractionation" annotation: short vertical arrow targeting the
  # Cerebellum-F bar (Female sits at the right side of the dodge group:
  # x = 3 + 0.85/4 ≈ 3.21). Label sits directly above the arrow tip.
  annotate("segment",
           x = 3.21, xend = 3.21,
           y = max(B_data$n_proteins) * 0.42,
           yend = qc_row$n_proteins[1] * 1.45,
           linewidth = 0.6, color = "grey25",
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("label",
           x = 3.21, y = max(B_data$n_proteins) * 0.50,
           label = "no fractionation",
           hjust = 0.5, size = 4.0, fontface = "bold", color = "grey20",
           fill = "white", label.size = 0.4) +
  scale_fill_manual(values = SEX_COLORS, na.value = "grey80") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.20))) +
  coord_cartesian(clip = "off") +
  labs(title = "Unique proteins per region * sex",
       x = NULL, y = "Unique proteins detected", fill = "Sex") +
  theme_poster() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.margin  = margin(10, 30, 10, 10))

# ----------------------------------------------------------------------------
# Panel C — 4-tissue Venn of QUANTIFIABLE proteins (proteotypic + doublet)
# Same definition used for the pyramid's "Quantifiable" tier.
# ----------------------------------------------------------------------------
prot_lists <- lapply(levels(peptides$tissue), function(tt) {
  peptides %>% filter(tissue == tt, proteotypic, doublet) %>%
    pull(protein_id) %>% unique()
})
names(prot_lists) <- levels(peptides$tissue)
cat("Quantifiable proteins per tissue:\n")
print(sapply(prot_lists, length))

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_grob <- venn.diagram(
  x = prot_lists[c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain")],
  category.names = c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain"),
  filename = NULL,
  fill = TISSUE_COLORS[c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain")],
  alpha = 0.45,
  lwd = 1.6,
  lty = "solid",
  col = TISSUE_COLORS[c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain")],
  cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.25,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = TISSUE_COLORS[c("Cerebrum", "BrainStem", "Cerebellum", "MockWholeBrain")],
  cat.dist = c(0.22, 0.22, 0.12, 0.12),
  margin = 0.10
)

n_shared4 <- length(Reduce(intersect, prot_lists))
panel_C <- cowplot::ggdraw() +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)) +
  cowplot::draw_label("Region overlap of quantifiable proteins",
                      x = 0.5, y = 0.97, size = 16, fontface = "bold",
                      colour = "black") +
  cowplot::draw_label(sprintf("%s shared across all 4 regions / proteotypic + doublet",
                              format(n_shared4, big.mark = ",")),
                      x = 0.5, y = 0.93, size = 11, colour = "grey30") +
  cowplot::draw_grob(grid::grobTree(venn_grob), x = 0, y = 0,
                    width = 1, height = 0.90)

# ----------------------------------------------------------------------------
# Panel D — Peptides per protein distribution by tissue (coverage depth)
# ----------------------------------------------------------------------------
D_data <- peptides %>%
  group_by(tissue, protein_id) %>%
  summarise(n_peptides = n_distinct(modified_peptide_sequence), .groups = "drop")

D_meds <- D_data %>%
  group_by(tissue) %>%
  summarise(med = median(n_peptides), n = n(), .groups = "drop")

panel_D <- ggplot(D_data, aes(x = tissue, y = n_peptides, fill = tissue)) +
  geom_violin(alpha = 0.85, scale = "width", color = "grey20", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.3, linewidth = 0.3) +
  geom_text(data = D_meds, aes(y = 1, label = sprintf("med=%d", med)),
            vjust = 1.6, size = 3.6, fontface = "bold") +
  geom_text(data = D_meds, aes(y = max(D_data$n_peptides),
                               label = sprintf("n=%s", format(n, big.mark = ","))),
            vjust = 1.0, size = 3.2, color = "grey25") +
  scale_fill_manual(values = TISSUE_COLORS, guide = "none") +
  scale_y_log10(labels = comma,
                breaks = c(1, 5, 10, 50, 100, 500),
                expand = expansion(mult = c(0.06, 0.05))) +
  labs(title = "Coverage depth: peptides per protein",
       x = NULL, y = "Peptides per protein (log10)") +
  theme_poster() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ----------------------------------------------------------------------------
# Combine panels (2 × 2)
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label(
    "Mouse Brain SysQuan - Resource Overview",
    fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62
  ) +
  draw_label(
    "TimsTOF DDA / 4 brain regions / 6 fractionated batches + 1 unfractionated (Cerebellum-F)",
    size = 13, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20
  ) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

top_row    <- plot_grid(panel_A, panel_B, ncol = 2,
                        labels = c("a", "b"), label_size = 18,
                        label_x = 0.01, label_y = 0.99,
                        rel_widths = c(1, 1))
bottom_row <- plot_grid(panel_C, panel_D, ncol = 2,
                        labels = c("c", "d"), label_size = 18,
                        label_x = 0.01, label_y = 0.99,
                        rel_widths = c(1.1, 1))

full <- plot_grid(header,
                  plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1.1)),
                  ncol = 1, rel_heights = c(0.07, 1))

save_poster(full, "FIG1_resource", width = 16, height = 12)

cat("\n=== FIGURE 1 written ===\n")
cat("  ../figures_poster/FIG1_resource.pdf\n")
cat("  ../figures_poster/FIG1_resource.png\n")
