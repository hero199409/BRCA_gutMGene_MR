# 01c_fig1_draft.R
# Figure 1 Mini Layout Draft: Parallel multi-cohort design & analytical boundaries
# Strictly avoids implying an unbroken "microbe -> metabolite -> gene -> therapy" causal proof.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Fig 1 Mini Layout Draft ===\n")

cards <- data.table(
  id   = c("prior", "mr_mib", "mr_met", "mr_cis", "scrna", "bulk", "ispy2", "aux_tx"),
  x0   = c(2,   24,  24,  24,  53.5, 53.5, 77,   77),
  x1   = c(20,  50,  50,  50,  74.5, 74.5, 98,   98),
  y0   = c(19,  47,  30,  11,  36,   11,   35,   11),
  y1   = c(58,  60,  44,  27,  60,   33,   60,   32),
  fill = c("#F8FAFC", "#EFF6FF", "#EFF6FF", "#EFF6FF", "#FEF2F2", "#F8FAFC", "#FEF2F2", "#F8FAFC"),
  bdr  = c("#94A3B8", "#3B82F6", "#3B82F6", "#3B82F6", "#DC2626", "#64748B", "#DC2626", "#94A3B8"),
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
    "Curated human microbe–\nmetabolite–gene links\n(61 candidate genes)\n\nSelects exposures &\ngenes for risk screens",
    "MiBioGen (64 taxa) vs\nBCAC (122,977 cases)\nReverse MR & DMP/FR02\nOutcome: BRCA risk",
    "47 metabolites vs BCAC\n30 taxa→metabolite\nStep-1 (min q = 0.429)\nOutcome: BRCA risk",
    "eQTLGen (40 cis GRCh37;\n36 harm. / 108 tests)\nColoc ±500 kb; 5 IL-8 pQTL\nOutcome: BRCA risk",
    "Wu (100k cells, 26 donors)\nE-MTAB-8107 (33k cells)\nGSE148673* (8.9k, no CXCL8)\n+ HPA IHC context",
    "5 cohorts (3,691 tumors)\nTCGA ER+ ΔC (n=770)\nTCGA TME (quanTIseq,\nMCP-counter) & GSEA",
    "Pretreatment FFPE\nAgilent 44K bulk RNA\nn=100 HER2− (29/71 vs 7/29)\nOutcome: pCR interaction",
    "GSE25066 / GSE32646 NACT\nMETABRIC / GEO Endocrine\nDepMap/GDSC2 (47×285;\nall FDR q >= 0.484)"
  )
)

banners <- data.table(
  x0 = c(2, 53.5, 77),
  x1 = c(50, 74.5, 98),
  y0 = c(62, 62, 62),
  y1 = c(66.5, 66.5, 66.5),
  fill = c("#1E3A8A", "#991B1B", "#991B1B"),
  label = c(
    "STREAM 1: Germline Genetic Screens (BRCA Risk)",
    "STREAM 2–3: Tumor Tissue",
    "STREAM 4: Therapy Response"
  )
)

arrows_dt <- data.table(
  x    = c(20.2, 20.2, 20.2),
  xend = c(23.7, 23.7, 23.7),
  y    = c(48,   37,   22),
  yend = c(53.5, 37,   19)
)

p_draft <- ggplot() +
  geom_rect(data = banners, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill), color = NA) +
  scale_fill_identity() +
  geom_text(data = banners, aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label),
            color = "white", size = 2.25, fontface = "bold") +
  geom_rect(data = cards, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = cards$fill, color = cards$bdr, linewidth = 0.65) +
  geom_text(data = cards, aes(x = (x0 + x1) / 2, y = y1 - 1.2, label = title),
            vjust = 1, size = 2.35, fontface = "bold", color = pal_rev$ink, lineheight = 1.02) +
  geom_text(data = cards, aes(x = (x0 + x1) / 2, y = y0 + (y1 - y0) * 0.36, label = body),
            vjust = 0.5, size = 2.00, color = pal_rev$subink, lineheight = 1.12) +
  geom_segment(data = arrows_dt, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
               linewidth = 0.65, color = pal_rev$pparg) +
  geom_segment(aes(x = 51.8, xend = 51.8, y = 10, yend = 66.5), linetype = "dashed", color = pal_rev$muted, linewidth = 0.5) +
  geom_segment(aes(x = 75.8, xend = 75.8, y = 10, yend = 66.5), linetype = "dashed", color = pal_rev$muted, linewidth = 0.5) +
  geom_rect(aes(xmin = 2, xmax = 98, ymin = 1.5, ymax = 8.5), fill = "#FFF7ED", color = "#D97706", linewidth = 0.55) +
  annotate("text", x = 50, y = 5.0,
           label = "Evidence & Measurement Boundary: Streams 1–4 analyze independent cohorts and endpoints (no joint serial-mediation chain is proven).\nI-SPY2 evaluates pretreatment bulk tumor RNA for pCR interaction; tumor myeloid localization of CXCL8 is established in independent scRNA-seq atlases.",
           size = 2.15, fontface = "bold", color = "#7C2D12", lineheight = 1.15) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 68), expand = FALSE) +
  labs(
    title = "Figure 1A (Layout Draft). Independent multi-cohort evidence streams, measurement layers, and endpoints",
    subtitle = "Arrows denote actual analytical flow (candidate prioritization within genetic risk screens);\nindependent tumor and trial streams are shown in parallel rather than as a proven causal chain."
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title    = element_text(face = "bold", size = 9.3, color = pal_rev$ink, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.3, color = pal_rev$subink, lineheight = 1.12, margin = margin(b = 4)),
    plot.margin   = margin(5, 6, 4, 6)
  )

save_rev_figure(p_draft, "Fig1_Layout_Draft", width_mm = 180, height_mm = 112, dpi_png = 350, out_dir = rev_prv_dir)
cat("Fig 1 Layout Draft polished complete.\n")
