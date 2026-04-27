# sysqun_mouse

Mouse Brain SysQuan poster analysis (CNPN 2026).

R pipeline that generates the seven poster-ready figures for the
mouse-brain proteome / phospho-proteome SysQuan dataset (4 brain regions,
2 sexes, TimsTOF DDA). Each figure is one self-contained script.

## Repository layout

```
sysqun_mouse/
├── analysis/
│   ├── 01_database_analysis.R       # data loader (sourced by 00_common.R)
│   ├── poster/
│   │   ├── 00_common.R              # shared loaders, palettes, helpers
│   │   ├── FIG1_resource.R          # resource overview
│   │   ├── FIG2_phospho.R           # phospho-proteome layer
│   │   ├── FIG3_tissue_sex.R        # tissue × sex multivariate
│   │   ├── FIG3b_sex_volcanoes.R    # proteome sex effect (volcanoes + GO)
│   │   ├── FIG3c_phospho_sex.R      # phospho-protein sex effect
│   │   ├── FIG3d_phosphopeptides.R  # phosphopeptide-level sex effect
│   │   └── FIG4_enrichment.R        # GO + driver-protein boxplots
│   ├── sysq_mouse_figs/             # final poster figures (PDF + PNG)
│   └── tables/                      # CSV outputs from each FIG script
├── database_statistics_*.txt        # DB summary stats (mouse / all-species)
├── .gitignore
└── README.md
```

## Required data (not committed)

The peptide bucket cache (`peptide_bucket_cache_20260423_115224.csv`,
~400 MB) is gitignored. Place it at the **repo root** before running:

```
sysqun_mouse/
├── peptide_bucket_cache_20260423_115224.csv   <-- here
└── analysis/
```

`analysis/01_database_analysis.R` auto-discovers any
`peptide_bucket_cache*.csv` in the parent directory and filters to the
four mouse-brain tissues.

## Running the pipeline

From `analysis/poster/`:

```bash
cd analysis/poster
Rscript FIG1_resource.R
Rscript FIG2_phospho.R
Rscript FIG3_tissue_sex.R
Rscript FIG3b_sex_volcanoes.R
Rscript FIG3c_phospho_sex.R
Rscript FIG3d_phosphopeptides.R
Rscript FIG4_enrichment.R
```

PDF + PNG land in `analysis/sysq_mouse_figs/`; CSV result tables land in
`analysis/tables/`.

## R dependencies

```r
install.packages(c(
  "tidyverse", "data.table", "cowplot", "patchwork", "ggrepel",
  "viridis", "scales", "RColorBrewer", "pheatmap",
  "VennDiagram", "ComplexUpset", "mixOmics", "BiocManager"
))
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))
```

## Methodology notes

- **Normalization.** Batch-median normalization on `log10(heavy_intensity)`
  (matches the human SysQuan paper convention). Stored as `log_heavy_norm`
  in `00_common.R`. Female batches have ~2× higher raw heavy intensity
  than males — a loading effect, not biology.
- **PLS-DA validation.** Every supervised PLS-DA reports a 999-permutation
  *p* (label shuffle vs observed between/within SS ratio). Subtitles flag
  non-significance.
- **Sample structure.** `n = 1` batch per (tissue × sex). No biological
  replicates. Cerebellum-Female was processed without fractionation
  (~20× fewer peptide rows) — included in coverage panels with a
  "no fractionation" annotation, excluded from multivariate analyses.
- **DE tests.** Peptide-level Welch within a protein, peptides as
  pseudo-replicates. P-values are descriptive ranking aids, not
  biological-replicate inference.
- **Phosphosite annotation.** `parse_phosphosite()` in `00_common.R`
  parses `[79.9663]` markers from modified peptide sequences and emits
  peptide-relative site labels (e.g. `Lrrc7 pS17`).
