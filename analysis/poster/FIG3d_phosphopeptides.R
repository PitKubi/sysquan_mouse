#!/usr/bin/env Rscript
# FIGURE 3d - Sex differential at PHOSPHOPEPTIDE resolution (CNPN 2026)
#
# Companion to FIG3c (which aggregates to phosphoprotein). Here each point
# is a single phosphopeptide. Per-peptide intensities are batch-median
# normalized via the `log_heavy_norm` column built in 00_common.R.
#
# Stat caveat: with n = 1 (peptide, sex) per region, a per-peptide Welch
# is undefined. We carry the protein-level Welch p-value down to each of
# its peptides for a stable y-axis. Multi-peptide proteins thus produce a
# vertical stack at one p, spread along x by per-peptide log2FC.
#
# Panels:
#   a. Cerebrum  Male vs Female phosphopeptide volcano
#   b. BrainStem Male vs Female phosphopeptide volcano
#   c. Per-peptide cross-region concordance (matched modified_peptide_sequence)
#   d. GO BP/MF/CC of phospho sex genes (Cer + BS + Cbl union)
#
# Run from analysis/poster/:  Rscript FIG3d_phosphopeptides.R

source("00_common.R")
suppressPackageStartupMessages({
  library(scales)
  library(ggrepel)
})

cat("\n=== FIGURE 3d - PHOSPHOPEPTIDE SEX VOLCANOES ===\n")

phos_pep <- peptides %>% filter(phospho)

# Protein-level Welch (used for the y-axis carrier)
sex_phos_cer_p <- sex_de_within_tissue(phos_pep, "Cerebrum",
                                       min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)
sex_phos_bs_p  <- sex_de_within_tissue(phos_pep, "BrainStem",
                                       min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)
sex_phos_cbl_p <- sex_de_within_tissue(phos_pep, "Cerebellum",
                                       min_n_per_sex = 1) %>%
  distinct(protein_id, .keep_all = TRUE)

# ----------------------------------------------------------------------------
# Build phosphopeptide-level M vs F log2FC + carry protein p-value
# ----------------------------------------------------------------------------
build_pep_volcano <- function(tissue_name, prot_de) {
  phos_pep %>%
    filter(tissue == tissue_name, !is.na(log_heavy_norm), !is.na(sex)) %>%
    group_by(protein_id, gene, modified_peptide_sequence, sex) %>%
    summarise(value = mean(log_heavy_norm, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sex, values_from = value) %>%
    filter(!is.na(Male), !is.na(Female)) %>%
    mutate(
      log2FC = (Male - Female) * log2(10),
      site_label = phosphosite_label(gene, modified_peptide_sequence)
    ) %>%
    left_join(prot_de %>% dplyr::select(protein_id,
                                        p_protein = p_value,
                                        adj_p_protein = adj_p),
              by = "protein_id") %>%
    mutate(
      direction = case_when(
        is.na(p_protein) | p_protein >= 0.05 ~ "n.s.",
        log2FC > 0 ~ "Up in Male",
        log2FC < 0 ~ "Up in Female",
        TRUE       ~ "n.s."
      ),
      direction = factor(direction, levels = c("Up in Male",
                                               "Up in Female", "n.s."))
    )
}

cer_pep <- build_pep_volcano("Cerebrum",  sex_phos_cer_p)
bs_pep  <- build_pep_volcano("BrainStem", sex_phos_bs_p)

cat(sprintf("Cerebrum  phosphopeptides: %d / with p: %d\n",
            nrow(cer_pep), sum(!is.na(cer_pep$p_protein))))
cat(sprintf("BrainStem phosphopeptides: %d / with p: %d\n",
            nrow(bs_pep),  sum(!is.na(bs_pep$p_protein))))

# ----------------------------------------------------------------------------
# Volcano helper - phosphopeptide-level
# ----------------------------------------------------------------------------
volcano_pep <- function(d, region_name, n_label = 14) {
  d_plot <- d %>% filter(!is.na(log2FC),
                         !is.na(p_protein), p_protein > 0)
  # Label top |log2FC| with explicit phosphosite (e.g. "Lrrc7 pS14").
  d_lab  <- d_plot %>%
    filter(!is.na(site_label), abs(log2FC) > 0.4) %>%
    arrange(desc(abs(log2FC))) %>%
    distinct(site_label, .keep_all = TRUE) %>%
    slice_head(n = n_label)
  ggplot(d_plot, aes(x = log2FC, y = -log10(p_protein), color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.55) +
    geom_point(data = d_lab, size = 2.3, alpha = 1) +
    geom_text_repel(data = d_lab, aes(label = site_label),
                    size = 3.3, color = "grey15", fontface = "italic",
                    box.padding = 0.45, point.padding = 0.3,
                    min.segment.length = 0.1, segment.color = "grey55",
                    max.overlaps = Inf, seed = 1) +
    scale_color_manual(values = DIRECTION_COLORS, drop = FALSE) +
    labs(title = sprintf("%s phosphopeptides: Male vs Female", region_name),
         subtitle = sprintf("%d phosphopeptides / x = peptide log2FC (norm.), y = protein-level -log10(p) / top |log2FC| labelled",
                            nrow(d_plot)),
         x = expression(log[2]~"(M/F), per phosphopeptide"),
         y = expression(-log[10]~"(protein p)"),
         color = NULL) +
    theme_poster() +
    theme(legend.position = "right")
}

panel_A <- volcano_pep(cer_pep, "Cerebrum")
panel_B <- volcano_pep(bs_pep,  "BrainStem")

# ----------------------------------------------------------------------------
# Panel C - cross-region peptide concordance (matched peptide sequence)
# ----------------------------------------------------------------------------
pep_match <- inner_join(
  cer_pep %>% dplyr::select(protein_id, gene, modified_peptide_sequence,
                            site_label,
                            log2FC_Cer = log2FC, p_Cer = p_protein),
  bs_pep  %>% dplyr::select(protein_id, modified_peptide_sequence,
                            log2FC_BS  = log2FC, p_BS  = p_protein),
  by = c("protein_id", "modified_peptide_sequence")
) %>%
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

cat(sprintf("Matched phosphopeptides Cer+BS: %d / concordant sig: %d\n",
            nrow(pep_match), sum(pep_match$concordant, na.rm = TRUE)))

# Phospho-sex concordance is rarely p-significant in this n=1 design, so we
# label the top per-peptide |mean log2FC| outliers regardless of class so
# the panel still tells a story.
lab_set <- pep_match %>%
  filter(!is.na(site_label)) %>%
  mutate(mean_abs = (abs(log2FC_Cer) + abs(log2FC_BS)) / 2) %>%
  arrange(desc(mean_abs)) %>%
  distinct(site_label, .keep_all = TRUE) %>%
  slice_head(n = 14)

class_pal <- c(
  "Shared - Up in Male"   = SEX_COLORS[["Male"]],
  "Shared - Up in Female" = SEX_COLORS[["Female"]],
  "Region-specific"       = "#fdae61",
  "n.s."                  = "grey80"
)

axis_lim <- max(abs(c(pep_match$log2FC_Cer, pep_match$log2FC_BS)),
                na.rm = TRUE)

panel_C <- ggplot(pep_match,
                  aes(x = log2FC_Cer, y = log2FC_BS, color = class)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.4) +
  geom_point(data = filter(pep_match, class == "n.s."),
             size = 1.0, alpha = 0.4) +
  geom_point(data = filter(pep_match, class != "n.s."),
             size = 2.2, alpha = 0.85) +
  geom_text_repel(data = lab_set, aes(label = site_label),
                  size = 3.3, color = "grey15", fontface = "italic",
                  box.padding = 0.45, point.padding = 0.3,
                  min.segment.length = 0.1, segment.color = "grey55",
                  max.overlaps = Inf, seed = 2) +
  scale_color_manual(values = class_pal, drop = FALSE) +
  scale_x_continuous(limits = c(-axis_lim, axis_lim)) +
  scale_y_continuous(limits = c(-axis_lim, axis_lim)) +
  coord_fixed() +
  labs(title = "Per-peptide cross-region concordance",
       subtitle = "Matched by modified peptide sequence; top |mean log2FC| labelled (n=1 per group, no biological replicates)",
       x = expression(log[2]~"(M/F) - Cerebrum (peptide)"),
       y = expression(log[2]~"(M/F) - BrainStem (peptide)"),
       color = NULL) +
  theme_poster() +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# Panel D - GO BP/MF/CC enrichment of phosphopeptide-derived sex genes
# Use top |log2FC| at the peptide level, mapped back to genes.
# ----------------------------------------------------------------------------
loose_pep_set <- function(d, top_n = 600) {
  d %>% filter(!is.na(log2FC), is.finite(log2FC),
               !is.na(gene), gene != "") %>%
    arrange(desc(abs(log2FC))) %>% slice_head(n = top_n)
}

pep_union <- bind_rows(
  loose_pep_set(cer_pep)  %>% mutate(region = "Cerebrum"),
  loose_pep_set(bs_pep)   %>% mutate(region = "BrainStem"),
  loose_pep_set(build_pep_volcano("Cerebellum", sex_phos_cbl_p))
    %>% mutate(region = "Cerebellum")
)

sex_genes_by_dir <- pep_union %>%
  group_by(gene) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(direction = ifelse(mean_log2FC > 0, "Up in Male", "Up in Female"))

cat(sprintf("Phosphopeptide-derived sex genes: %d up-M / %d up-F\n",
            sum(sex_genes_by_dir$direction == "Up in Male"),
            sum(sex_genes_by_dir$direction == "Up in Female")))

phos_universe <- phos_pep %>%
  filter(!is.na(gene), gene != "") %>% pull(gene) %>% unique()

go_pep <- bind_rows(lapply(c("Up in Male", "Up in Female"), function(dir) {
  syms <- sex_genes_by_dir %>% filter(direction == dir) %>%
    pull(gene) %>% unique()
  enrich_gene_set_all(syms, label = dir, fdr_cap = 0.25, top_n = 5,
                      universe = phos_universe)
}))

panel_D <- plot_go_3panels(
  go_pep,
  title    = "GO enrichment - phosphopeptide-derived sex genes",
  subtitle = sprintf("BP / MF / CC top 5 each; %d up-M, %d up-F genes (Cer+BS+Cbl, top 600 |log2FC| per region)",
                     sum(sex_genes_by_dir$direction == "Up in Male"),
                     sum(sex_genes_by_dir$direction == "Up in Female"))
)

# ----------------------------------------------------------------------------
# Combine
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Sex Effect at Phosphopeptide Resolution",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("Per-peptide volcanoes + concordance + GO / batch-median normalized heavy / n=1 batch per (tissue,sex) - no biological replicates",
             size = 12, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20) +
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

save_poster(full, "FIG3d_phosphopeptides", width = 17, height = 17)

write_csv(pep_match %>% arrange(desc(concordant), p_Cer + p_BS),
          file.path(TAB_DIR, "FIG3d_phosphopeptide_concordance.csv"))
write_csv(go_pep, file.path(TAB_DIR, "FIG3d_GO_phosphopeptides.csv"))

cat("\n=== FIGURE 3d written ===\n")
cat(sprintf("  %s/FIG3d_phosphopeptides.pdf\n", FIG_DIR))
cat(sprintf("  %s/FIG3d_phosphopeptides.png\n", FIG_DIR))
