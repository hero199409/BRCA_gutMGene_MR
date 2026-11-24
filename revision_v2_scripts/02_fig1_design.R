# 02_fig1_design.R
# Figure 1 (Final Clean Minimalist Version): Independent multi-cohort evidence streams,
# measurement layers, and endpoints. Zero prose subtitle/footer inside graphic.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Final Fig 1 (Study Design & Evidence Boundaries) ===\n")

cards <- data.table(
  id   = c("prior", "mr_mib", "mr_met", "mr_cis", "scrna", "bulk", "ispy2", "aux_tx"),
  x0   = c(1.5, 23.5, 23.5, 23.5, 53.0, 53.0, 76.8, 76.8),
  x1   = c(19.5, 49.5, 49.5, 49.5, 74.2, 74.2, 98.5, 98.5),
  y0   = c(19,  47,   30,   11,   36,   11,   35,   11),
  y1   = c(58,  60,   44,   27,   60,   33,   60,   32),
  fill = c("#F8FAFC", "#EFF6FF", "#EFF6FF", "#EFF6FF", "#FEF2F2", "#F8FAFC", "#FEF2F2", "#F8FAFC"),
  bdr  = c("#94A3B8", "#2B5B84", "#2B5B84", "#2B5B84", "#C83E2B", "#64748B", "#C83E2B", "#94A3B8"),
  title = c(
    "Experimental Prior\ngutMGene v1.0 / v2.0",
    "1a. Microbiome MR",
    "1b. Metabolite MR",
    "1c. Regulatory cis-QTL MR",
    "2. Tumor scRNA-seq Atlases\n(Cell Localization)",
    "3. Bulk Cohorts & TME\n(Prognosis & Pathways)",
    "4a. Primary Trial Test\nI-SPY2 (GSE173839)",
    "4b. Auxiliary Checks\nNACT, Endo & GDSC2"
  ),
  body = c(
    "Curated human microbe\u2013\nmetabolite\u2013gene links\n(61 candidate genes)\n\nPrioritizes exposures &\ngenes for risk screens",
    "MiBioGen (64 taxa) vs\nBCAC (122,977 cases)\nReverse MR & DMP/FR02\nOutcome: BRCA risk",
    "47 metabolites vs BCAC\n30 taxa\u2192metabolite\nStep-1 (min q = 0.429)\nOutcome: BRCA risk",
    "eQTLGen (40 cis GRCh37;\n36 harm. / 108 tests)\nColoc \u00b1500 kb; 5 IL-8 pQTL\nOutcome: BRCA risk",
    "Wu (100k cells, 26 donors)\nE-MTAB-8107 (33k cells)\nGSE148673* (8.9k, no CXCL8)\n+ HPA IHC context",
    "5 cohorts (3,691 tumors)\nTCGA ER+ \u0394C (n = 770)\nTCGA TME (quanTIseq,\nMCP-counter) & GSEA",
    "Pretreatment FFPE\nAgilent 44K bulk RNA\nn=100 HER2\u2212 (29/71 vs 7/29)\nOutcome: pCR interaction",
    "GSE25066 / GSE32646 NACT\nMETABRIC / GEO Endocrine\nDepMap/GDSC2 (47\u00d7285;\nall FDR q \u2265 0.484)"
  )
)

banners <- data.table(
  x0 = c(1.5, 53.0, 76.8),
  x1 = c(49.5, 74.2, 98.5),
  y0 = c(62, 62, 62),
  y1 = c(66.8, 66.8, 66.8),
  fill = c("#1E3A8A", "#991B1B", "#991B1B"),
  label = c(
    "STREAM 1: Germline Genetic Screens (BRCA Risk)",
    "STREAM 2\u20133: Tumor Tissue",
    "STREAM 4: Therapy Response"
  )
)

arrows_dt <- data.table(
  x    = c(19.8, 19.8, 19.8),
  xend = c(23.2, 23.2, 23.2),
  y    = c(48,   37,   22),
  yend = c(53.5, 37,   19)
)

p_fig1 <- ggplot() +
  geom_rect(data = banners, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill), color = NA) +
  scale_fill_identity() +
  geom_text(data = banners, aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label),
            color = "white", size = 2.35, fontface = "bold") +
  geom_rect(data = cards, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = cards$fill, color = cards$bdr, linewidth = 0.65) +
  geom_text(data = cards, aes(x = (x0 + x1) / 2, y = y1 - 1.2, label = title),
            vjust = 1, size = 2.40, fontface = "bold", color = pal_rev$ink, lineheight = 1.03) +
  geom_text(data = cards, aes(x = (x0 + x1) / 2, y = y0 + (y1 - y0) * 0.36, label = body),
            vjust = 0.5, size = 2.05, color = pal_rev$subink, lineheight = 1.13) +
  geom_segment(data = arrows_dt, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
               linewidth = 0.65, color = pal_rev$pparg) +
  geom_segment(aes(x = 51.2, xend = 51.2, y = 10, yend = 66.8), linetype = "dashed", color = pal_rev$muted, linewidth = 0.5) +
  geom_segment(aes(x = 75.5, xend = 75.5, y = 10, yend = 66.8), linetype = "dashed", color = pal_rev$muted, linewidth = 0.5) +
  geom_rect(aes(xmin = 1.5, xmax = 98.5, ymin = 1.5, ymax = 8.4), fill = "#FFF7ED", color = "#D97706", linewidth = 0.55) +
  annotate("text", x = 50, y = 4.95,
           label = "Measurement & Inference Boundary: Streams 1\u20134 evaluate independent cohorts and endpoints (no joint causal chain is tested).\nI-SPY2 profiles pretreatment bulk tumor RNA; myeloid CXCL8 localization is established in independent scRNA-seq atlases.",
           size = 2.20, fontface = "bold", color = "#7C2D12", lineheight = 1.15) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 68.5), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(3, 3, 3, 3))

save_rev_figure(p_fig1, "Fig1_Study_Design_and_Evidence_Boundaries", width_mm = 180, height_mm = 106, dpi_png = 350)
cat("Fig 1 complete.\n")
