# 03_supp_figsS1toS4.R
# Supplementary Figures S1-S4 in minimalist top-journal style (180 mm width, 350 dpi PNG + vector PDF).
# S1: Upstream MR sensitivity scatter & leave-one-out
# S2: CXCL8 OS heterogeneity LOO & same-rank cross-cohort microbiome MR meta-analysis
# S3: External scRNA UMAPs (E-MTAB-8107, GSE148673), HPA single-cell nCPM, GSE9195 ER+ DMFS
# S4: Relocated auxiliary clinical therapy (Endocrine, NACT pCR/DRFS) & GDSC2 pharmacogenomics

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

# ==============================================================================
# FIG S1: Upstream MR Sensitivity & Leave-One-Out
# ==============================================================================
cat("=== Building Fig S1 (Upstream MR Sensitivity & LOO) ===\n")

single_dt <- fread(file.path(res_dir, "planA_mr_single_snp_wald.csv"))
sens_dt   <- fread(file.path(res_dir, "planA_mr_sensitivity_summary.csv"))
loo_dt    <- fread(file.path(res_dir, "planA_mr_leaveoneout.csv"))

target_s1 <- list(
  list(eid = "ebi-a-GCST90016960", oid = "ieu-a-1127", tag = "A1", ltag = "B1", lbl = "Adlercreutzia \u2192 ER+"),
  list(eid = "ebi-a-GCST90016987", oid = "ieu-a-1128", tag = "A2", ltag = "B2", lbl = "Desulfovibrio \u2192 ER\u2212"),
  list(eid = "ebi-a-GCST90017039", oid = "ieu-a-1126", tag = "A3", ltag = "B3", lbl = "Parabacteroides \u2192 Overall"),
  list(eid = "ebi-a-GCST90016951", oid = "ieu-a-1127", tag = "A4", ltag = "B4", lbl = "Ruminococcaceae \u2192 ER+"),
  list(eid = "ebi-a-GCST90016989", oid = "ieu-a-1128", tag = "A5", ltag = "B5", lbl = "Dorea \u2192 ER\u2212"),
  list(eid = "met-d-Acetate",      oid = "ieu-a-1126", tag = "A6", ltag = "B6", lbl = "Acetate \u2192 Overall")
)

mr_cols <- c(
  "Inverse variance weighted" = pal_rev$cxcl8,
  "Weighted median"           = pal_rev$teal,
  "MR Egger"                  = pal_rev$pparg,
  "Weighted mode"             = pal_rev$amber
)

sc_plots  <- list()
loo_plots <- list()
for (tk in target_s1) {
  panel_label <- sub(" \u2192 ", "\n\u2192 ", tk$lbl, fixed = TRUE)
  sd <- copy(single_dt[id.exposure == tk$eid & id.outcome == tk$oid])
  flip <- sd$beta.exposure < 0
  sd[, bx := ifelse(flip, -beta.exposure, beta.exposure)]
  sd[, by := ifelse(flip, -beta.outcome, beta.outcome)]
  ms <- copy(sens_dt[id.exposure == tk$eid & id.outcome == tk$oid & method %in% names(mr_cols)])
  ms[, intercept := ifelse(method == "MR Egger" & !is.na(egger_intercept), egger_intercept, 0)]
  x_rng <- range(sd$bx, na.rm = TRUE)
  ms[, `:=`(x1 = x_rng[1] * 0.85, x2 = x_rng[2] * 1.06,
            y1 = intercept + b * x_rng[1] * 0.85, y2 = intercept + b * x_rng[2] * 1.06)]
  ms[, method := factor(method, levels = names(mr_cols))]

  p_sc <- ggplot(sd, aes(x = bx, y = by)) +
    geom_hline(yintercept = 0, color = pal_rev$border, linetype = "dotted") +
    geom_errorbar(aes(ymin = by - 1.96 * se.outcome, ymax = by + 1.96 * se.outcome),
                  color = pal_rev$muted, width = 0, linewidth = 0.35) +
    geom_errorbar(aes(xmin = bx - 1.96 * se.exposure, xmax = bx + 1.96 * se.exposure),
                  orientation = "y", color = pal_rev$muted, width = 0, linewidth = 0.35) +
    geom_segment(data = ms, aes(x = x1, xend = x2, y = y1, yend = y2, color = method), linewidth = 0.82) +
    geom_point(size = 1.75, color = pal_rev$ink) +
    scale_color_manual(values = mr_cols, breaks = names(mr_cols), name = NULL, drop = FALSE) +
    labs(title = sprintf("%s  %s", tk$tag, panel_label), x = "SNP effect (exposure)", y = "SNP effect (BRCA)") +
    theme_rev(8.0) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 8.1, lineheight = 1.02, margin = margin(b = 3)),
          axis.text = element_text(size = 6.9),
          axis.title = element_text(size = 7.4))
  sc_plots[[tk$lbl]] <- p_sc

  ld <- copy(loo_dt[id.exposure == tk$eid & id.outcome == tk$oid])
  ld[, is_all := (SNP == "All")]
  ld <- ld[order(is_all, -b)]
  ld[, snp_ord := factor(SNP, levels = unique(rev(SNP)))]
  ld[, y_num := as.numeric(snp_ord)]
  z_lo <- ld[y_num %% 2 == 1]

  p_lo <- ggplot(ld, aes(x = b, y = snp_ord, color = is_all)) +
    geom_rect(data = z_lo, aes(ymin = y_num - 0.46, ymax = y_num + 0.46, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, fill = pal_rev$zebra, color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
    geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), orientation = "y", width = 0.24, linewidth = 0.5) +
    geom_point(size = 1.75) +
    scale_color_manual(values = c("FALSE" = pal_rev$slate, "TRUE" = pal_rev$cxcl8), guide = "none") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
    labs(title = sprintf("%s  %s", tk$ltag, panel_label), x = expression("IVW " * beta), y = NULL) +
    theme_rev(8.0) +
    theme(plot.title = element_text(size = 8.1, lineheight = 1.02, margin = margin(b = 3)),
          axis.text.y = element_text(size = 6.4),
          axis.text.x = element_text(size = 6.9),
          axis.title.x = element_text(size = 7.4))
  loo_plots[[tk$lbl]] <- p_lo
}

s1_leg   <- cowplot::get_legend(sc_plots[[1]] + theme(legend.position = "bottom",
                                                      legend.text = element_text(size = 8.0),
                                                      legend.key.width = unit(5.5, "mm")))
sc_clean <- lapply(sc_plots, function(p) p + theme(legend.position = "none"))
s1_grid  <- wrap_plots(c(sc_clean, loo_plots), ncol = 3)
fig_s1   <- s1_grid / patchwork::wrap_elements(full = s1_leg) + plot_layout(heights = c(1.0, 0.045))
save_rev_figure(fig_s1, "FigS1_Upstream_MR_Sensitivity", width_mm = 180, height_mm = 198, dpi_png = 350)
rm(fig_s1, sc_plots, sc_clean, loo_plots, s1_grid); gc()

# ==============================================================================
# FIG S2: CXCL8 OS Heterogeneity LOO & Same-Rank Cross-Cohort Microbiome Meta
# ==============================================================================
cat("=== Building Fig S2 (Heterogeneity LOO & Microbiome Meta) ===\n")

het_map <- c(
  "All_5_Cohorts_Overall_OS"        = "All 5 cohorts (OS)",
  "LOO_Excluding_METABRIC"          = "Excl. METABRIC",
  "LOO_Excluding_GSE42568"          = "Excl. GSE42568",
  "LOO_Excluding_TCGA_BRCA"         = "Excl. TCGA-BRCA",
  "LOO_Excluding_GSE20685"          = "Excl. GSE20685",
  "LOO_Excluding_GSE7390"           = "Excl. GSE7390",
  "Microarray_4Cohorts_Overall_OS"  = "Microarray cohorts only (k = 4)",
  "GEO_3Cohorts_Affy_OS"            = "Affymetrix GEO (k = 3)",
  "GEO_Large_GSE20685_GSE7390_OS"   = "Large Affymetrix GEO (GSE20685 + GSE7390)",
  "ER_neg_Stratum_TCGA_METABRIC_OS" = "TCGA + METABRIC ER\u2212 stratum"
)

het_dt <- fread(file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))
het_dt[, `:=`(
  row_label = ifelse(analysis %in% names(het_map), unname(het_map[analysis]), gsub("_", " ", analysis)),
  effect    = RE_HR, low = RE_LCI, high = RE_UCI,
  stats     = sprintf("%.2f [%.2f, %.2f] | %.3g | %.0f%%", RE_HR, RE_LCI, RE_UCI, RE_pval, I2_pct),
  col_group = fifelse(analysis == "All_5_Cohorts_Overall_OS", pal_rev$cxcl8, pal_rev$pparg)
)]
ps2a <- table_forest_rev(
  het_dt,
  title_txt = "A   CXCL8 overall survival heterogeneity and leave-one-out",
  x_lab = "Random-effects HR (log scale)", null_val = 1,
  stats_header = "HR [95% CI] | P | I\u00b2",
  widths = c(1.68, 1.08, 1.16), log_x = TRUE, base_size = 8.2
)

cc_meta <- fread(file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))[analysis_scope == "same_rank"]
stopifnot(nrow(cc_meta) > 0L, !anyDuplicated(cc_meta, by = c("core_taxon", "outcome_label")))
cc_meta[, outcome_short := clean_subtype_rev(outcome_label)]
cc_meta[, outcome_order := match(outcome_short, c("Overall", "ER+", "ER\u2212"))]
setorder(cc_meta, core_taxon, outcome_order)
cc_meta[, `:=`(
  row_label = sprintf("%s | %s | k=%d", core_taxon, outcome_short, k_cohorts),
  effect    = RE_OR, low = RE_LCI, high = RE_UCI,
  stats     = sprintf("%.3f [%.3f, %.3f] | %.3g | %.3f | %.0f%%", RE_OR, RE_LCI, RE_UCI, RE_pval, RE_FDR, I2_pct),
  col_group = fifelse(RE_pval < 0.05, pal_rev$pparg, pal_rev$slate)
)]
ps2b <- table_forest_rev(
  cc_meta,
  title_txt = "B   Cross-cohort microbiome MR meta-analysis (same rank)",
  x_lab = "Random-effects OR (log scale)", null_val = 1,
  stats_header = "OR [95% CI] | P | BH q | I\u00b2",
  widths = c(1.38, 0.98, 1.56), log_x = TRUE, base_size = 8.2
)

fig_s2 <- ps2a / ps2b + plot_layout(heights = c(1.0, 1.62))
save_rev_figure(fig_s2, "FigS2_HeterogeneityLOO_and_MicrobiomeMeta", width_mm = 180, height_mm = 184, dpi_png = 350)
rm(fig_s2, ps2a, ps2b); gc()

# ==============================================================================
# FIG S3: External scRNA UMAPs, HPA Single-Cell nCPM, & GSE9195 ER+ DMFS
# ==============================================================================
cat("=== Building Fig S3 (External scRNA + HPA nCPM + GSE9195 DMFS) ===\n")

sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))

plot_umap_ct_rev <- function(df, pal, title_txt) {
  d <- copy(df)[!is.na(UMAP_1) & !is.na(UMAP_2)]
  d[, Lineage := clean_lineage_rev(Celltype_major)]
  ct_lvls <- d[, .N, by = Lineage][order(-N)]$Lineage
  d[, Lineage := factor(Lineage, levels = ct_lvls)]
  meds <- d[, .(x = median(UMAP_1), y = median(UMAP_2), n = .N), by = Lineage][n >= max(60, 0.01 * nrow(d))]
  xr <- range(d$UMAP_1); yr <- range(d$UMAP_2)
  x0 <- xr[1]; y0 <- yr[1]; dx <- diff(xr) * 0.16; dy <- diff(yr) * 0.16

  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = Lineage)) +
    geom_point(size = 0.16, alpha = 0.65, show.legend = FALSE) +
    ggrepel::geom_label_repel(
      data = meds, aes(x = x, y = y, label = Lineage),
      inherit.aes = FALSE, size = 2.0, fontface = "bold", color = pal_rev$ink,
      fill = alpha("white", 0.82), label.padding = unit(0.08, "lines"),
      linewidth = 0.2, segment.size = 0.22, segment.color = pal_rev$muted,
      box.padding = 0.15, max.overlaps = 20, seed = 42
    ) +
    annotate("segment", x = x0, xend = x0 + dx, y = y0, yend = y0,
             arrow = arrow(length = unit(1.1, "mm"), type = "closed"), linewidth = 0.35, color = pal_rev$ink) +
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + dy,
             arrow = arrow(length = unit(1.1, "mm"), type = "closed"), linewidth = 0.35, color = pal_rev$ink) +
    annotate("text", x = x0 + dx * 0.5, y = y0 - dy * 0.25, label = "UMAP1", size = 1.85, color = pal_rev$slate) +
    annotate("text", x = x0 - dx * 0.25, y = y0 + dy * 0.5, label = "UMAP2", size = 1.85, angle = 90, color = pal_rev$slate) +
    scale_color_manual(values = pal) +
    labs(title = title_txt) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(size = 8.4, face = "bold", color = pal_rev$ink, margin = margin(b = 2)),
          plot.margin = margin(2, 2, 2, 2))
}

plot_umap_gene_rev <- function(df, gene_col, title_txt, pal_hi) {
  d <- copy(df)[!is.na(UMAP_1) & !is.na(UMAP_2)]
  d[, val := get(gene_col)]
  setorder(d, val)
  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = val)) +
    geom_point(size = 0.16, alpha = 0.80) +
    scale_color_gradientn(colors = c("#E2E8F0", "#CBD5E1", pal_hi), name = NULL,
                          guide = guide_colorbar(barwidth = unit(2.2, "mm"), barheight = unit(13, "mm"))) +
    labs(title = title_txt) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(size = 8.4, face = "bold", color = pal_rev$ink, margin = margin(b = 2)),
          legend.position = "right",
          legend.text = element_text(size = 6.5, color = pal_rev$slate),
          legend.margin = margin(0, 0, 0, -4),
          plot.margin = margin(2, 2, 2, 2))
}

em_df <- copy(sc_cache[["BRCA_EMTAB8107"]])
ct_em <- clean_lineage_rev(em_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_em <- setNames(colorRampPalette(c(pal_rev$pparg, pal_rev$cxcl8, pal_rev$teal, "#0284C7",
                                      pal_rev$coral, pal_rev$slate, pal_rev$amber, "#7C3AED", "#EC4899"))(length(ct_em)), ct_em)
ps3a1 <- plot_umap_ct_rev(em_df, pal_em, "A1  E-MTAB-8107 (n=33,043)")
ps3a2 <- plot_umap_gene_rev(em_df, "PPARG", "A2  PPARG", pal_rev$pparg)
ps3a3 <- plot_umap_gene_rev(em_df, "CXCL8", "A3  CXCL8", pal_rev$cxcl8)

g14_df <- copy(sc_cache[["BRCA_GSE148673"]])
ct_g14 <- clean_lineage_rev(g14_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_g14 <- setNames(colorRampPalette(c(pal_rev$pparg, pal_rev$cxcl8, pal_rev$teal, "#0284C7",
                                       pal_rev$coral, pal_rev$slate, pal_rev$amber, "#7C3AED"))(length(ct_g14)), ct_g14)
ps3b1 <- plot_umap_ct_rev(g14_df, pal_g14, "B1  GSE148673 (n=10,359)")
ps3b2 <- plot_umap_gene_rev(g14_df, "PPARG", "B2  PPARG", pal_rev$pparg)
ps3b3 <- plot_umap_gene_rev(g14_df, "IL1B",  "B3  IL1B",  pal_rev$cxcl8)

sc_hpa <- fread(file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
sc_keep <- c("Neutrophils", "Macrophages", "Fibroblasts", "Adipocytes",
             "Lymphatic endothelial cells", "Vascular endothelial cells",
             "Breast secretory cells", "T-cells", "B-cells")
sc_sub <- sc_hpa[gene_symbol %in% c("CXCL8", "PPARG") & cell_type %in% sc_keep]
sc_sub[, cell_clean := tools::toTitleCase(gsub("-", " ", cell_type))]
sc_sub[, log10_ncpm := log10(nCPM + 1)]
sc_ord <- sc_sub[, .(mx = max(log10_ncpm)), by = cell_clean][order(mx)]$cell_clean
sc_sub[, cell_clean := factor(cell_clean, levels = sc_ord)]
sc_sub[, y_num := as.numeric(cell_clean)]
z_s3c <- data.table(y_num = seq_along(sc_ord))[y_num %% 2 == 1]

ps3c <- ggplot(sc_sub, aes(x = log10_ncpm, y = cell_clean, fill = gene_symbol)) +
  geom_rect(data = z_s3c, aes(ymin = y_num - 0.46, ymax = y_num + 0.46, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = pal_rev$zebra, color = NA) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66) +
  geom_text(aes(label = ifelse(nCPM > 0, format(round(nCPM, 0), big.mark = ","), "0"),
                x = log10_ncpm + 0.06),
            position = position_dodge(width = 0.74), hjust = 0, size = 2.05, color = pal_rev$slate) +
  scale_fill_manual(values = c("CXCL8" = pal_rev$cxcl8, "PPARG" = pal_rev$pparg), name = NULL) +
  coord_cartesian(xlim = c(0, 4.45)) +
  labs(title = "C   HPA single-cell nCPM", x = expression(bold(log[10](nCPM + 1))), y = NULL) +
  theme_rev(8.2) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        axis.text.y = element_text(size = 7.2))

th_pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))
g91_km <- copy(th_pc$gse9195)[!is.na(t_dmfs) & t_dmfs > 0 & !is.na(e_dmfs)]
g91_km[, Group := ifelse(PPARG >= median(PPARG), "High", "Low")]
km_s3d <- plot_km_rev(
  g91_km, "t_dmfs", "e_dmfs", "Group",
  "D   GSE9195 ER+ DMFS (PPARG)", "Years", "DMFS probability",
  "PPARG", c("Low", "High"), c(pal_rev$cxcl8, pal_rev$pparg),
  break_by = 3, max_x = 12, base_size = 8.0
)

figs3_r1 <- (ps3a1 | ps3a2 | ps3a3)
figs3_r2 <- (ps3b1 | ps3b2 | ps3b3)
figs3_r3 <- patchwork::wrap_elements(full = (ps3c | km_s3d) + plot_layout(widths = c(1.06, 1.04)))
fig_s3   <- figs3_r1 / figs3_r2 / figs3_r3 + plot_layout(heights = c(0.95, 0.95, 1.18))
save_rev_figure(fig_s3, "FigS3_scRNA_ExternalCohorts_and_GSE9195", width_mm = 180, height_mm = 190, dpi_png = 350)
rm(fig_s3, figs3_r1, figs3_r2, figs3_r3, ps3a1, ps3a2, ps3a3, ps3b1, ps3b2, ps3b3, ps3c, km_s3d, sc_cache, em_df, g14_df); gc()

# ==============================================================================
# FIG S4: Relocated Auxiliary Clinical Therapy & Subtype-Adjusted GDSC2
# ==============================================================================
cat("=== Building Fig S4 (Auxiliary Clinical Therapy & GDSC2) ===\n")

th_dt <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))

# S4A: PPARG & ER+ endocrine therapy
endo_dt <- copy(th_dt[gene == "PPARG" & cohort %in% c("METABRIC", "GSE9195", "GSE17705") & stratum == "ER_pos"])
endo_dt[, arm_clean := fcase(
  treatment_arm == "Hormone_Therapy_YES", "Endocrine (+)",
  treatment_arm == "Hormone_Therapy_NO",  "Endocrine (\u2212)",
  treatment_arm == "Tamoxifen_5yr",       "Tamoxifen",
  default = treatment_arm
)]
endo_dt[, `:=`(
  row_label = sprintf("%s | %s | %s (n=%s)", cohort, arm_clean, endpoint, format(n_samples, big.mark = ",", trim = TRUE)),
  effect    = estimate, low = lci95, high = uci95,
  stats     = sprintf("%.2f [%.2f, %.2f] | %.3g", estimate, lci95, uci95, pval),
  col_group = fifelse(pval < 0.05, pal_rev$pparg, pal_rev$slate)
)]
ps4a <- table_forest_rev(
  endo_dt, "A   PPARG and ER+ endocrine therapy outcomes", "HR (log scale)", 1,
  "HR [95% CI] | Cox P", c(1.95, 1.35, 1.08), log_x = TRUE, base_size = 8.2
)

# S4B: CXCL8 & NACT pCR (single-line labels to prevent Old Fig 8B collision)
nact_pcr_dt <- copy(th_dt[gene == "CXCL8" & cohort %in% c("GSE25066", "GSE32646") &
                          endpoint == "pCR_vs_RD" & stratum %in% c("Overall", "ER_neg", "TNBC")])
nact_pcr_dt[, `:=`(
  row_label = sprintf("%s %s (n=%d, %d pCR)", cohort, clean_subtype_rev(stratum), n_samples, n_events),
  effect    = estimate, low = lci95, high = uci95,
  stats     = sprintf("%.2f | %.3g", auc, pval),
  col_group = fifelse(pval < 0.05, pal_rev$cxcl8, pal_rev$slate)
)]
ps4b <- table_forest_rev(
  nact_pcr_dt, "B   CXCL8 and NACT pCR", "OR (log scale)", 1,
  "AUC | P", c(1.62, 1.02, 0.56), log_x = TRUE, x_breaks = c(0.6, 1, 2), base_size = 8.0
)

# S4C: GSE25066 NACT DRFS KM
g25_km <- copy(th_pc$gse25066)[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs)]
g25_km[, Group := ifelse(CXCL8 >= median(CXCL8), "High", "Low")]
km_s4c <- plot_km_rev(
  g25_km, "t_drfs", "e_drfs", "Group",
  "C   GSE25066 NACT DRFS (CXCL8)", "Years", "DRFS probability",
  "CXCL8", c("Low", "High"), c(pal_rev$pparg, pal_rev$cxcl8),
  break_by = 2, max_x = 8, base_size = 7.9
)

# S4D: GDSC2 subtype-adjusted drug response (CXCL8 & PPARG side-by-side)
sub_d <- fread(file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
adj   <- sub_d[gene_symbol %in% c("CXCL8", "PPARG") & is.finite(subtype_adjusted_pval)]
stopifnot(!anyDuplicated(adj, by = c("gene_symbol", "DRUG_ID")))
bh_chk <- adj[, .(ok = isTRUE(all.equal(subtype_adjusted_FDR, p.adjust(subtype_adjusted_pval, method = "BH"), tolerance = 1e-10)) &&
                       all(subtype_adjusted_family_n == .N)), by = gene_symbol]
stopifnot(nrow(bh_chk) == 2L, all(bh_chk$ok))

make_gdsc_forest <- function(gene, example_id, x_range, title_txt, col_pt) {
  gd <- adj[gene_symbol == gene]
  key <- gd[is_key_drug == TRUE & DRUG_ID != example_id][order(subtype_adjusted_pval, DRUG_ID)]
  ex  <- gd[DRUG_ID == example_id]
  d <- rbind(key[seq_len(5L)], ex)
  setorder(d, subtype_adjusted_pval, DRUG_ID)
  d[, row_label := gsub("_", "-", DRUG_NAME, fixed = TRUE)]
  d[duplicated(row_label) | duplicated(row_label, fromLast = TRUE),
    row_label := sprintf("%s (%s)", row_label, DRUG_ID)]
  d[, `:=`(effect = subtype_adjusted_beta, low = subtype_adjusted_ci_low, high = subtype_adjusted_ci_high,
           stats = sprintf("%.3f", subtype_adjusted_FDR))]
  table_forest_rev(d, title_txt, "Adjusted \u03b2 (95% CI)", 0, "BH q",
                   c(1.15, 1.38, 0.42), log_x = FALSE, x_range = x_range,
                   point_col = col_pt, base_size = 8.0)
}

adj_range <- adj[is_key_drug == TRUE, range(c(subtype_adjusted_ci_low, subtype_adjusted_ci_high), finite = TRUE)]
ps4d1 <- make_gdsc_forest("CXCL8", 1079L, adj_range, "D1  GDSC2 CXCL8 (m=285)", pal_rev$cxcl8)
ps4d2 <- make_gdsc_forest("PPARG", 1200L, adj_range, "D2  GDSC2 PPARG (m=285)", pal_rev$pparg)

# S4E1 & S4E2: GDSC2 scatter with boxed non-overlapping labels + single shared legend
gb_mat <- readRDS(file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))
gb_mat[, Subtype := factor(Subtype_Group, levels = c("ER_pos", "HER2_pos", "TNBC"), labels = c("ER+", "HER2+", "TNBC"))]

spearman_box <- function(gene_name, drug_id, d, target_label) {
  rho <- sub_d[gene_symbol == gene_name & DRUG_ID == drug_id]
  d <- d[is.finite(get(gene_name)) & is.finite(LN_IC50)]
  stopifnot(nrow(rho) == 1L, nrow(d) == rho$n_all, nrow(d[Subtype_Group == rho$target_subtype]) == rho$n_target_subtype)
  sprintf("All: \u03c1 = %.3f (P = %.3f)\n%s: \u03c1 = %.3f (P = %.3f)",
          rho$all_spearman_rho, rho$all_spearman_pval,
          target_label, rho$within_subtype_rho, rho$within_subtype_pval)
}

sub_cols <- c("ER+" = pal_rev$pparg, "HER2+" = pal_rev$amber, "TNBC" = pal_rev$cxcl8)

d_e1 <- gb_mat[DRUG_ID == 1079L]
e1_lbl <- spearman_box("CXCL8", 1079L, d_e1, "TNBC")
ps4e1 <- ggplot(d_e1, aes(x = CXCL8, y = LN_IC50, color = Subtype)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = pal_rev$ink, linetype = "dashed", linewidth = 0.65, alpha = 0.14) +
  geom_smooth(data = d_e1[Subtype_Group == "TNBC"], method = "lm", se = FALSE, color = pal_rev$cxcl8, linewidth = 0.80) +
  geom_point(size = 2.0, alpha = 0.88) +
  annotate("label", x = max(d_e1$CXCL8, na.rm = TRUE), y = max(d_e1$LN_IC50, na.rm = TRUE) + 0.55,
           label = e1_lbl, hjust = 1, vjust = 1, size = 2.25, color = pal_rev$ink,
           fill = alpha("white", 0.90), linewidth = 0.25, label.padding = unit(0.18, "lines")) +
  scale_color_manual(values = sub_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  labs(title = "E1  CXCL8 vs dasatinib",
       x = expression(bold(CXCL8~(log[2]~TPM + 1))), y = expression(bold(Dasatinib~ln~IC[50]))) +
  theme_rev(8.2) + theme(legend.position = "bottom")

d_e2 <- gb_mat[DRUG_ID == 1200L]
e2_lbl <- spearman_box("PPARG", 1200L, d_e2, "ER+")
ps4e2 <- ggplot(d_e2, aes(x = PPARG, y = LN_IC50, color = Subtype)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = pal_rev$ink, linetype = "dashed", linewidth = 0.65, alpha = 0.14) +
  geom_smooth(data = d_e2[Subtype_Group == "ER_pos"], method = "lm", se = FALSE, color = pal_rev$pparg, linewidth = 0.80) +
  geom_point(size = 2.0, alpha = 0.88) +
  annotate("label", x = min(d_e2$PPARG, na.rm = TRUE), y = max(d_e2$LN_IC50, na.rm = TRUE) + 0.45,
           label = e2_lbl, hjust = 0, vjust = 1, size = 2.25, color = pal_rev$ink,
           fill = alpha("white", 0.90), linewidth = 0.25, label.padding = unit(0.18, "lines")) +
  scale_color_manual(values = sub_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  labs(title = "E2  PPARG vs fulvestrant",
       x = expression(bold(PPARG~(log[2]~TPM + 1))), y = expression(bold(Fulvestrant~ln~IC[50]))) +
  theme_rev(8.2) + theme(legend.position = "bottom")

s4_leg <- cowplot::get_legend(ps4e1 + theme(legend.position = "bottom", legend.text = element_text(size = 8.0)))
row_s4_2 <- patchwork::wrap_elements(full = (ps4b | km_s4c) + plot_layout(widths = c(1.16, 1.0)))
row_s4_3 <- patchwork::wrap_elements(full = (ps4d1 | ps4d2))
row_s4_4 <- patchwork::wrap_elements(full = ((ps4e1 + theme(legend.position = "none")) |
                                             (ps4e2 + theme(legend.position = "none"))))

fig_s4 <- ps4a / row_s4_2 / row_s4_3 / row_s4_4 / patchwork::wrap_elements(full = s4_leg) +
  plot_layout(heights = c(1.05, 1.38, 1.25, 1.18, 0.16))
save_rev_figure(fig_s4, "FigS4_Auxiliary_Therapy_and_GDSC2", width_mm = 180, height_mm = 236, dpi_png = 350)

cat("All Supplementary Figures S1-S4 complete.\n")
