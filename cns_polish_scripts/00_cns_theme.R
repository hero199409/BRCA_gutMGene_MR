# 00_cns_theme.R
# Unified CNS / Nature-journal visual system for BRCA gutMGene MR Manuscript
# Isolated directory: 06_figures_cns_polish
# NEVER modifies anything inside 05_revision_v2.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(ggrepel)
  library(ggbeeswarm)
  library(ggrastr)
  library(ragg)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(jpeg)
  library(png)
  library(magick)
  library(grid)
  library(scales)
})

# ------------------------------------------------------------------------------
# Read-only data paths (strictly read-only)
# ------------------------------------------------------------------------------
root_dir  <- Sys.getenv("BRCA_PROJECT_ROOT", unset = ".")
cache_dir <- file.path(root_dir, "00_data/mr_cache")
bulk_dir  <- file.path(root_dir, "00_data/bulk_cohorts")
geo_dir   <- file.path(root_dir, "00_data/geo_cache")
sc_dir    <- file.path(root_dir, "00_data/scRNA_cache")
res_dir   <- file.path(root_dir, "02_results")
ihc_dir   <- file.path(root_dir, "03_figures/HPA_IHC_Images")

# ------------------------------------------------------------------------------
# Isolated output paths (only write here)
# ------------------------------------------------------------------------------
cns_dir     <- file.path(root_dir, "06_figures_cns_polish")
cns_scr_dir <- file.path(cns_dir, "01_scripts")
cns_fig_dir <- file.path(cns_dir, "03_figures")
dir.create(cns_scr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cns_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# NPG (Nature Publishing Group) refined soft palette
# ------------------------------------------------------------------------------
pal_cns <- list(
  # Primary signal (CXCL8 / positive)
  cxcl8       = "#E64B35",
  cxcl8_dark  = "#B83220",
  cxcl8_soft  = "#FDE8E5",
  # Control / negative (PPARG)
  pparg       = "#3C5488",
  pparg_light = "#4DBBD5",
  pparg_dark  = "#243456",
  pparg_soft  = "#E6F4F8",
  # Auxiliary highlight (immunotherapy / interaction / optimism-corrected)
  teal        = "#00A087",
  teal_dark   = "#007563",
  teal_soft   = "#E6F5F2",
  # Secondary accent
  salmon      = "#F39B7F",
  amber       = "#D97706",
  purple      = "#8491B4",
  # Neutral / structural tones
  ink         = "#111827",
  subink      = "#374151",
  slate       = "#6B7280",
  muted       = "#9CA3AF",
  ref_line    = "#CBD5E1",
  grid_line   = "#EEF2F6",
  zebra       = "#F8FAFC",
  border_soft = "#E2E8F0"
)

# ------------------------------------------------------------------------------
# Typography & theme system (Arial sans-serif, ultra-thin 0.35-0.40 pt axes,
# white frameless facet strips, ggtext markdown support for italic gene symbols)
# ------------------------------------------------------------------------------
cns_font <- "Arial"

theme_cns <- function(base_size = 8.0, base_family = cns_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Ultra-thin crisp axes (0.35-0.40 pt)
      axis.line         = element_line(color = pal_cns$ink, linewidth = 0.38),
      axis.ticks        = element_line(color = pal_cns$subink, linewidth = 0.35),
      axis.ticks.length = unit(2.2, "pt"),
      axis.title.x      = ggtext::element_markdown(size = base_size - 0.2, color = pal_cns$ink,
                                                   margin = margin(t = 4)),
      axis.title.y      = ggtext::element_markdown(size = base_size - 0.2, color = pal_cns$ink,
                                                   margin = margin(r = 4)),
      axis.text.x       = ggtext::element_markdown(size = base_size - 0.8, color = pal_cns$subink),
      axis.text.y       = ggtext::element_markdown(size = base_size - 0.8, color = pal_cns$subink),
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      # Frameless pure-white facet strips with clean typography
      strip.background  = element_blank(),
      strip.text        = ggtext::element_markdown(face = "bold", size = base_size - 0.4,
                                                   color = pal_cns$ink, hjust = 0,
                                                   margin = margin(t = 3, b = 3, l = 0, r = 2)),
      strip.placement   = "outside",
      # Clear hierarchy: bold large tag + plain smaller subtitle-style panel title
      plot.tag          = element_text(face = "bold", size = base_size + 3.5, color = pal_cns$ink),
      plot.tag.position = c(0, 1),
      plot.title        = ggtext::element_markdown(face = "plain", size = base_size + 0.6,
                                                   color = pal_cns$ink, hjust = 0,
                                                   margin = margin(t = 0.5, b = 5, l = 13)),
      plot.title.position = "plot",
      plot.subtitle     = ggtext::element_markdown(size = base_size - 0.8, color = pal_cns$slate,
                                                   hjust = 0, margin = margin(b = 3, l = 13)),
      plot.caption      = ggtext::element_markdown(size = base_size - 1.1, color = pal_cns$slate,
                                                   hjust = 0, margin = margin(t = 3)),
      plot.caption.position = "plot",
      # Borderless transparent legends
      legend.title      = ggtext::element_markdown(face = "bold", size = base_size - 0.7, color = pal_cns$ink),
      legend.text       = ggtext::element_markdown(size = base_size - 0.8, color = pal_cns$subink),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key        = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.margin       = margin(5, 7, 5, 5)
    )
}

theme_cns_grid <- function(base_size = 8.0, base_family = cns_font) {
  theme_cns(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major.y = element_line(color = pal_cns$grid_line, linewidth = 0.32),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# ------------------------------------------------------------------------------
# Save helper (PNG 350 DPI via ragg + vector PDF via cairo_pdf)
# ------------------------------------------------------------------------------
save_cns_figure <- function(p, name, width_mm = 180, height_mm = 180, dpi_png = 350, out_dir = cns_fig_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  w_in <- width_mm / 25.4
  h_in <- height_mm / 25.4
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = w_in, height = h_in, dpi = dpi_png, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = w_in, height = h_in, device = cairo_pdf, bg = "white")
  cat(sprintf("Saved %-48s | %.0fx%.0f mm | PNG=%.2f MB | PDF=%.2f MB\n",
              name, width_mm, height_mm,
              file.info(png_path)$size / 1048576,
              file.info(pdf_path)$size / 1048576))
}

# ------------------------------------------------------------------------------
# Shared formatting helpers (strictly preserving statistical values)
# ------------------------------------------------------------------------------
format_p_cns <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-16 && p == 0) return("< 2.2e-16")
  if (p < 0.001) return(formatC(p, format = "e", digits = 1))
  sprintf("%.3f", p)
}

fmt_minus <- function(x) {
  gsub("-", "\u2212", x, fixed = TRUE)
}

clean_subtype_cns <- function(x) {
  x <- gsub("BRCA_", "", x)
  x <- gsub("ER_pos", "ER+", x)
  x <- gsub("ER_neg", "ER\u2212", x)
  x <- gsub("HER2_pos", "HER2+", x)
  x
}

clean_lineage_cns <- function(x) {
  m <- c(
    "Mono/Macro"     = "Myeloid",
    "T_cells"        = "T cells",
    "CD4Tconv"       = "CD4+ T",
    "CD8T"           = "CD8+ T",
    "CD8Tex"         = "Exh. CD8+ T",
    "Tprolif"        = "Prolif. T",
    "B"              = "B cells",
    "Plasma"         = "Plasma",
    "Epithelial"     = "Norm. epithelial",
    "Malignant"      = "Malignant",
    "Myofibroblasts" = "Myofibroblasts",
    "Fibroblasts"    = "Fibroblasts",
    "Endothelial"    = "Endothelial",
    "PVL"            = "PVL",
    "Mast"           = "Mast",
    "DC"             = "DC"
  )
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}
