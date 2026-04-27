#!/usr/bin/env Rscript
# FIGURE 4 — Functional enrichment (CNPN 2026)
# Panels:
#   a. GO BP enrichment per region (one-vs-rest, full proteome)
#   b. GO BP enrichment by sex direction (concordant Cerebrum+BrainStem)
#   c. GO BP enrichment of phospho region markers (Cerebrum/BrainStem/Cerebellum)
#   d. Driver-gene boxplots (5-6 genes that anchor the pathway story)
#
# Run from analysis/poster/:  Rscript FIG4_enrichment.R

source("00_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(scales)
})

cat("\n=== FIGURE 4 - FUNCTIONAL ENRICHMENT ===\n")

# ----------------------------------------------------------------------------
# Common: build a "universe" of detected genes (with a known mouse symbol)
# ----------------------------------------------------------------------------
universe_syms <- peptides %>%
  filter(!is.na(gene), gene != "") %>%
  pull(gene) %>% unique()
cat(sprintf("Universe genes: %d\n", length(universe_syms)))

# ----------------------------------------------------------------------------
# Helper: GO BP enrichment from a gene-symbol vector vs background universe
# ----------------------------------------------------------------------------
enrich_gene_set <- function(gene_syms, label, fdr_cap = 0.10, top_n = 4) {
  gene_syms <- unique(gene_syms[!is.na(gene_syms) & gene_syms != ""])
  if (length(gene_syms) < 3) {
    return(tibble(label = label, Description = NA, p.adjust = NA,
                  GeneRatio = NA, Count = NA))
  }
  ego <- tryCatch(
    enrichGO(gene = gene_syms,
             universe = universe_syms,
             OrgDb = org.Mm.eg.db,
             keyType = "SYMBOL",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = fdr_cap,
             readable = FALSE,
             minGSSize = 5, maxGSSize = 1000),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(tibble(label = label, Description = NA, p.adjust = NA,
                  GeneRatio = NA, Count = NA))
  }
  res <- as.data.frame(ego) %>%
    mutate(label = label,
           gene_ratio_n = as.numeric(sub("/.*", "", GeneRatio)) /
                          as.numeric(sub(".*/", "", GeneRatio))) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    dplyr::select(label, Description, p.adjust, gene_ratio_n, Count)
  res
}

# ----------------------------------------------------------------------------
# Panel A — Per-region GO BP enrichment (full proteome)
# ----------------------------------------------------------------------------
TISSUES_ALL <- levels(peptides$tissue)
tissue_de_all <- bind_rows(lapply(TISSUES_ALL, function(tt) {
  tissue_de_one_vs_rest(peptides, tt, channel = "heavy_intensity", min_n = 3)
})) %>% distinct(tissue, protein_id, .keep_all = TRUE)

# Use rank-based top-N (top 200 by log2FC, requiring p<0.05) — standard approach
# when biological replication is limited (peptide-level p-values are
# pseudo-replicate-driven, so signed log2FC is the more honest ranker).
go_tissue <- bind_rows(lapply(TISSUES_ALL, function(tt) {
  syms <- tissue_de_all %>%
    filter(tissue == tt, !is.na(p_value), p_value < 0.05, log2FC > 0,
           !is.na(gene), gene != "") %>%
    arrange(desc(log2FC)) %>%
    slice_head(n = 200) %>%
    pull(gene) %>% unique()
  cat(sprintf("  Tissue %-15s: %d up-regulated genes (top 200 by log2FC, p<0.05)\n",
              tt, length(syms)))
  enrich_gene_set(syms, label = tt, fdr_cap = 0.10, top_n = 4)
}))
go_tissue <- go_tissue %>% filter(!is.na(Description))

panel_A <- ggplot(go_tissue,
                  aes(x = factor(label, levels = TISSUES_ALL),
                      y = stringr::str_wrap(Description, 40),
                      size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradientn(
    colors = c("#1f77b4", "#2ca02c", "#d62728"),
    name = expression(-log[10]~"adj.p")) +
  scale_size_continuous(range = c(2.5, 8), name = "Gene\noverlap") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "GO Biological Process - per-region (proteome)",
       subtitle = "Top 4 BP terms per region (one-vs-rest, FDR<0.10)",
       x = NULL, y = NULL) +
  theme_poster() +
  theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "right")

# ----------------------------------------------------------------------------
# Panel B — Sex GO BP enrichment (Cerebrum + BrainStem concordant)
# ----------------------------------------------------------------------------
sex_de_cer <- sex_de_within_tissue(peptides, "Cerebrum",  channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)
sex_de_bs  <- sex_de_within_tissue(peptides, "BrainStem", channel = "heavy_intensity") %>%
  distinct(protein_id, .keep_all = TRUE)

sex_concord <- inner_join(
  sex_de_cer %>% dplyr::select(protein_id, gene, log2FC_Cer = log2FC,
                               adj_p_Cer = adj_p, p_Cer = p_value),
  sex_de_bs  %>% dplyr::select(protein_id, log2FC_BS = log2FC,
                               adj_p_BS = adj_p,  p_BS = p_value),
  by = "protein_id"
) %>%
  filter(!is.na(gene), gene != "",
         sign(log2FC_Cer) == sign(log2FC_BS),
         pmin(p_Cer, p_BS, na.rm = TRUE) < 0.05) %>%   # p (peptide-level), not FDR
  mutate(
    mean_log2FC = (log2FC_Cer + log2FC_BS) / 2,
    direction   = ifelse(mean_log2FC > 0, "Up in Male", "Up in Female")
  )

cat(sprintf("Sex-concordant gene set sizes: Up-M=%d, Up-F=%d\n",
            sum(sex_concord$direction == "Up in Male"),
            sum(sex_concord$direction == "Up in Female")))

go_sex <- bind_rows(lapply(c("Up in Male", "Up in Female"), function(dir) {
  syms <- sex_concord %>% filter(direction == dir) %>% pull(gene) %>% unique()
  enrich_gene_set(syms, label = dir, fdr_cap = 0.25, top_n = 4)  # loosen for tiny lists
}))
go_sex <- go_sex %>% filter(!is.na(Description))

panel_B <- if (nrow(go_sex) > 0) {
  ggplot(go_sex,
         aes(x = factor(label, levels = c("Up in Male", "Up in Female")),
             y = stringr::str_wrap(Description, 40),
             size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradientn(colors = c("#1f77b4", "#2ca02c", "#d62728"),
                          name = expression(-log[10]~"adj.p")) +
    scale_size_continuous(range = c(2.5, 8), name = "Gene\noverlap") +
    labs(title = "GO Biological Process - sex direction (proteome)",
         subtitle = sprintf("Concordant Cerebrum + BrainStem, %d up-M / %d up-F genes",
                            sum(sex_concord$direction == "Up in Male"),
                            sum(sex_concord$direction == "Up in Female")),
         x = NULL, y = NULL) +
    theme_poster() +
    theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
          legend.position = "right")
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, size = 6, color = "grey40",
                      label = "No GO BP terms passed FDR<0.25\n(sex gene lists too small)") +
    theme_void(base_size = 14) +
    labs(title = "GO Biological Process - sex direction") +
    theme(plot.title = element_text(face = "bold", size = 14))
}

# ----------------------------------------------------------------------------
# Panel C — Phospho region marker GO BP (Cer/BS/Cbl)
# ----------------------------------------------------------------------------
PHOS_TISSUES <- c("Cerebrum", "BrainStem", "Cerebellum")
phos_pep <- peptides %>% filter(tissue %in% PHOS_TISSUES, phospho)
phos_tissue_de <- bind_rows(lapply(PHOS_TISSUES, function(tt) {
  tissue_de_one_vs_rest(phos_pep, tt, channel = "heavy_intensity", min_n = 3)
})) %>% distinct(tissue, protein_id, .keep_all = TRUE)

# Phospho one-vs-rest p-values are sparse (small phospho subset), so take
# top 100 by log2FC ignoring p-value, restrict to log2FC > 0.5.
go_phos <- bind_rows(lapply(PHOS_TISSUES, function(tt) {
  syms <- phos_tissue_de %>%
    filter(tissue == tt, !is.na(log2FC), log2FC > 0.5,
           !is.na(gene), gene != "") %>%
    arrange(desc(log2FC)) %>%
    slice_head(n = 100) %>%
    pull(gene) %>% unique()
  cat(sprintf("  Phospho-marker %-12s: %d up-regulated genes (top 100 by log2FC, log2FC>0.5)\n",
              tt, length(syms)))
  enrich_gene_set(syms, label = tt, fdr_cap = 0.25, top_n = 4)
}))
go_phos <- go_phos %>% filter(!is.na(Description))

panel_C <- if (nrow(go_phos) > 0) {
  ggplot(go_phos,
         aes(x = factor(label, levels = PHOS_TISSUES),
             y = stringr::str_wrap(Description, 40),
             size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradientn(colors = c("#1f77b4", "#2ca02c", "#d62728"),
                          name = expression(-log[10]~"adj.p")) +
    scale_size_continuous(range = c(2.5, 8), name = "Gene\noverlap") +
    labs(title = "GO Biological Process - phospho-region markers",
         subtitle = "Cerebrum, BrainStem, Cerebellum (FDR<0.25; small input lists)",
         x = NULL, y = NULL) +
    theme_poster() +
    theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
          legend.position = "right")
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, size = 6, color = "grey40",
                      label = "No phospho-marker GO BP terms passed FDR<0.25") +
    theme_void(base_size = 14) +
    labs(title = "GO Biological Process - phospho regions") +
    theme(plot.title = element_text(face = "bold", size = 14))
}

# ----------------------------------------------------------------------------
# Panel D — Driver-gene boxplots
# Choose 6 driver genes: top tissue marker per region + top sex driver
# ----------------------------------------------------------------------------
top_per_region <- tissue_de_all %>%
  filter(!is.na(adj_p), log2FC > 0, !is.na(gene), gene != "") %>%
  group_by(tissue) %>% arrange(adj_p, desc(log2FC), .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup()

top_sex <- sex_concord %>%
  arrange(desc(abs(log2FC_Cer + log2FC_BS) / 2)) %>%
  slice_head(n = 2)

driver_genes <- unique(c(top_per_region$gene, top_sex$gene))
cat("Driver genes:", paste(driver_genes, collapse = ", "), "\n")

driver_data <- peptides %>%
  filter(gene %in% driver_genes, !is.na(heavy_intensity), heavy_intensity > 0,
         !batch %in% MV_BAD_BATCHES, !is.na(sex)) %>%
  mutate(log_int = log10(heavy_intensity)) %>%
  group_by(gene, tissue, sex, batch) %>%
  summarise(mean_log10 = mean(log_int, na.rm = TRUE), .groups = "drop")

# Order facets from low -> high overall protein expression (median across all batches)
gene_order <- driver_data %>%
  group_by(gene) %>%
  summarise(med = median(mean_log10, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>% pull(gene)
driver_data$gene <- factor(driver_data$gene, levels = gene_order)

panel_D <- ggplot(driver_data, aes(x = tissue, y = mean_log10, fill = sex)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8),
               width = 0.7, color = "grey25", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.8), size = 2.0,
             shape = 21, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = SEX_COLORS, na.value = "grey80") +
  facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  labs(title = "Driver protein expression across regions and sex",
       subtitle = "Top region markers + top sex-concordant proteins / facets sorted low -> high median expression",
       x = NULL, y = expression(log[10]~"(heavy intensity)")) +
  theme_poster() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        strip.text  = element_text(face = "italic", size = 12),
        plot.subtitle = element_text(size = 10, color = "grey30"))

# ----------------------------------------------------------------------------
# Combine: 2x2 grid (A | B) / (C | D)
# ----------------------------------------------------------------------------
header <- ggdraw() +
  draw_label("Mouse Brain SysQuan - Functional Enrichment",
             fontface = "bold", size = 22, hjust = 0.5, x = 0.5, y = 0.62) +
  draw_label("GO Biological Process enrichment of regional, sex, and phospho drivers",
             size = 13, colour = "grey30", hjust = 0.5, x = 0.5, y = 0.20) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

row1 <- plot_grid(panel_A, panel_B, ncol = 2,
                  labels = c("a", "b"), label_size = 18,
                  label_x = 0.01, label_y = 0.99,
                  rel_widths = c(1.1, 1))
row2 <- plot_grid(panel_C, panel_D, ncol = 2,
                  labels = c("c", "d"), label_size = 18,
                  label_x = 0.01, label_y = 0.99,
                  rel_widths = c(1, 1.05))

full <- plot_grid(header,
                  plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1)),
                  ncol = 1, rel_heights = c(0.06, 1))

save_poster(full, "FIG4_enrichment", width = 17, height = 13)

write_csv(go_tissue, file.path(TAB_DIR, "FIG4_GO_tissue.csv"))
write_csv(go_sex,    file.path(TAB_DIR, "FIG4_GO_sex.csv"))
write_csv(go_phos,   file.path(TAB_DIR, "FIG4_GO_phospho.csv"))

cat("\n=== FIGURE 4 written ===\n")
cat("  ../figures_poster/FIG4_enrichment.pdf\n")
cat("  ../figures_poster/FIG4_enrichment.png\n")
