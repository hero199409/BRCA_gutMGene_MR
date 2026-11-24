# 06_fig6_pparg_discrim_cns.R
# Figure 6 (Nature Medicine / Cancer Cell high-impact composite upgrade):
#   Panel A: UpSet-style Covariate Inclusion Matrix x Multivariable Forest Plot + HR/P table
#   Panel B: 5-tier Model Discrimination Ladder (C-index +/- SE) with delta-AIC & LRT P badges
#   Panel C: Dual-layer Time-dependent IPCW AUC trajectories + delta-AUC diverging lollipop
#   Panel D: Free-scale Bootstrap delta-C Optimism Shrinkage Vector plot with 95% CI bands
#
# Strictly preserves 100% of statistical values, sample sizes, CIs, and data sources from 05_revision_v2.
# Writes exclusively to 06_figures_cns_polish/03_figures.

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "06_figures_cns_polish/01_scripts/00_cns_theme.R"))

cat("=== Generating High-Impact CNS Fig 6 (4-Panel Composite, Polished) ===\n")

ad_dt  <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
paired <- fread(file.path(res_dir, "planA_erpos_delta_c_bootstrap.csv"))
roc_dt <- fread(file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))

coh_pal   <- c("TCGA-BRCA" = pal_cns$pparg, "METABRIC" = pal_cns$cxcl8)
coh_short <- c(TCGA_ERpos = "TCGA-BRCA", METABRIC_ERpos = "METABRIC")
coh_strip <- c(
  TCGA_ERpos     = "<span style='color:#3C5488;'>TCGA-BRCA</span> <span style='font-weight:normal;color:#6B7280;'>(n = 770, 94 deaths)</span>",
  METABRIC_ERpos = "<span style='color:#E64B35;'>METABRIC</span> <span style='font-weight:normal;color:#6B7280;'>(n = 1,127, 626 deaths)</span>"
)
comp_labels <- c(Clinical = "Clinical", Clinical_Adipose = "Clinical + adipose")

fmt_panel_hdr <- function(letter, text) {
  sprintf("<span style='font-size:11.5pt;font-weight:bold;color:#111827;'>%s</span>\u2002\u2002<span style='font-size:8.6pt;color:#111827;'>%s</span>",
          letter, text)
}

# ==============================================================================
# PANEL A: UpSet-Style Covariate Inclusion Matrix x Forest Plot + HR & P Table
# ==============================================================================
ad_sub <- ad_dt[endpoint == "OS" & cohort %chin% c("TCGA_ERpos", "METABRIC_ERpos")]
ad_sub[, cohort_s := factor(coh_short[cohort], levels = c("TCGA-BRCA", "METABRIC"))]
ad_sub[, row_in_coh := seq_len(.N), by = cohort]
ad_sub[, y_pos := ifelse(cohort == "TCGA_ERpos", 12.2 - row_in_coh, 6.0 - row_in_coh)]
ad_sub[, sig := pval < 0.05]
ad_sub[, model_tag := fcase(
  model == "M1_Univariate", "M1",
  startsWith(model, "M2_"), "M2",
  model == "M3_Adj_Adipose+Stromal", "M3",
  model == "M3_Adj_Age+Stage", "M3",
  model == "M4_Adj_Age+Stage", "M4",
  model == "M4_Adj_Age+Stage+Grade", "M4",
  startsWith(model, "M5_"), "M5"
)]
ad_sub[, hr_txt := sprintf("%.2f (%.2f\u2013%.2f)", HR, HR_LCI95, HR_UCI95)]
ad_sub[, p_txt  := ifelse(pval < 0.001, formatC(pval, format = "e", digits = 1), sprintf("%.3f", pval))]
ad_sub[, p_col  := ifelse(sig, ifelse(cohort == "TCGA_ERpos", pal_cns$pparg, pal_cns$cxcl8), pal_cns$slate)]
ad_sub[, pt_col := ifelse(cohort == "TCGA_ERpos", pal_cns$pparg, pal_cns$cxcl8)]
ad_sub[, is_zebra := (row_in_coh %% 2L == 1L)]

cov_cols   <- c("PPARG", "Adipose", "Stromal", "Age", "Stage", "Grade")
cov_labels <- c("*PPARG*", "*ADIPOQ*/*PLIN1*", "Stromal", "Age", "Stage", "Grade")

ad_sub[, `:=`(
  c_PPARG   = TRUE,
  c_Adipose = grepl("Adipose|Full", model),
  c_Stromal = grepl("Stromal", model),
  c_Age     = grepl("Age|Full", model),
  c_Stage   = grepl("Stage|Full", model),
  c_Grade   = grepl("Grade", model)
)]

mat_dt <- melt(
  ad_sub[, .(cohort, cohort_s, model_tag, y_pos, pt_col,
             PPARG = c_PPARG, Adipose = c_Adipose, Stromal = c_Stromal,
             Age = c_Age, Stage = c_Stage, Grade = c_Grade)],
  id.vars = c("cohort", "cohort_s", "model_tag", "y_pos", "pt_col"),
  variable.name = "cov", value.name = "included"
)
mat_dt[, cov_idx := match(as.character(cov), cov_cols)]

seg_dt <- mat_dt[included == TRUE, .(x_min = min(cov_idx), x_max = max(cov_idx),
                                     pt_col = pt_col[1]), by = .(y_pos)]

coh_hdr_a <- data.table(
  y_pos = c(12.02, 5.82),
  lbl = c(
    "<span style='color:#3C5488;'><b>TCGA-BRCA</b></span> <span style='color:#6B7280;'>(n = 770, 94 deaths)</span>",
    "<span style='color:#E64B35;'><b>METABRIC</b></span> <span style='color:#6B7280;'>(n = 1,127, 626 deaths)</span>"
  )
)

y_lim_a <- c(0.45, 12.60)

# A1: Left UpSet Covariate Inclusion Matrix (clip = "off" so cohort headers never clip)
p6a_mat <- ggplot(mat_dt, aes(x = cov_idx, y = y_pos)) +
  geom_rect(data = ad_sub[is_zebra == TRUE],
            aes(ymin = y_pos - 0.44, ymax = y_pos + 0.44),
            xmin = 0.25, xmax = 6.55, fill = pal_cns$zebra, inherit.aes = FALSE) +
  geom_segment(data = seg_dt[x_max > x_min],
               aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos, color = pt_col),
               linewidth = 0.70, alpha = 0.55, inherit.aes = FALSE) +
  geom_point(data = mat_dt[included == FALSE],
             color = "#D8DEE9", fill = "#F1F5F9", shape = 21, size = 2.0, stroke = 0.40) +
  geom_point(data = mat_dt[included == TRUE],
             aes(fill = pt_col, color = pt_col), shape = 21, size = 2.25, stroke = 0.35) +
  geom_text(data = ad_sub, aes(x = 0.20, y = y_pos, label = model_tag),
            hjust = 1, size = 2.05, fontface = "bold", color = pal_cns$slate,
            family = cns_font, inherit.aes = FALSE) +
  ggtext::geom_richtext(data = coh_hdr_a, aes(x = -0.32, y = y_pos, label = lbl),
                        hjust = 0, vjust = 0.5, size = 2.15, fill = NA, label.color = NA,
                        label.padding = unit(0, "pt"), family = cns_font, inherit.aes = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = 1:6, labels = cov_labels, expand = c(0, 0), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.35, 6.65), ylim = y_lim_a, clip = "off") +
  labs(x = "Model covariates", y = NULL) +
  theme_cns(7.8) +
  theme(
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.line.x  = element_line(color = pal_cns$border_soft, linewidth = 0.35),
    axis.ticks.x = element_blank(),
    axis.text.x  = ggtext::element_markdown(angle = 35, hjust = 1, vjust = 1,
                                            size = 6.7, color = pal_cns$subink),
    axis.title.x = ggtext::element_markdown(size = 7.3, color = pal_cns$ink, margin = margin(t = 2)),
    plot.margin  = margin(2, 2, 2, 4)
  )

# A2: Middle Forest Plot (segmented vertical null line so it never crosses top text)
p6a_for <- ggplot(ad_sub, aes(x = HR, y = y_pos)) +
  annotate("rect", xmin = 0.54, xmax = 1.0, ymin = 0.55, ymax = 5.45,
           fill = "#F4F8FB", alpha = 0.85) +
  annotate("rect", xmin = 0.54, xmax = 1.0, ymin = 6.75, ymax = 11.65,
           fill = "#F4F8FB", alpha = 0.85) +
  geom_rect(data = ad_sub[is_zebra == TRUE],
            aes(ymin = y_pos - 0.44, ymax = y_pos + 0.44),
            xmin = 0.54, xmax = 1.28, fill = pal_cns$zebra, alpha = 0.55, inherit.aes = FALSE) +
  annotate("segment", x = 1, xend = 1, y = 0.55, yend = 5.45,
           linetype = "dashed", color = pal_cns$muted, linewidth = 0.42) +
  annotate("segment", x = 1, xend = 1, y = 6.75, yend = 11.65,
           linetype = "dashed", color = pal_cns$muted, linewidth = 0.42) +
  geom_linerange(aes(xmin = HR_LCI95, xmax = HR_UCI95, color = pt_col), linewidth = 0.60) +
  geom_point(aes(color = pt_col, fill = ifelse(sig, pt_col, "white")),
             shape = 21, size = 2.35, stroke = 0.80) +
  annotate("text", x = 0.96, y = 12.02, label = "\u2190 Lower risk", hjust = 1,
           size = 1.92, color = pal_cns$pparg, fontface = "italic", family = cns_font) +
  annotate("text", x = 1.04, y = 12.02, label = "Higher \u2192", hjust = 0,
           size = 1.92, color = pal_cns$slate, fontface = "italic", family = cns_font) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.2), labels = c("0.6", "0.8", "1.0", "1.2")) +
  scale_y_continuous(limits = y_lim_a, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.54, 1.28), clip = "off") +
  labs(x = "HR per native unit (95% CI)", y = NULL) +
  theme_cns(7.8) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = ggtext::element_markdown(size = 7.3, color = pal_cns$ink, margin = margin(t = 5)),
    plot.margin  = margin(2, 3, 2, 2)
  )

# A3: Right HR (95% CI) & P-value Table
p6a_tbl <- ggplot(ad_sub, aes(y = y_pos)) +
  geom_rect(data = ad_sub[is_zebra == TRUE],
            aes(ymin = y_pos - 0.44, ymax = y_pos + 0.44),
            xmin = 0, xmax = 1.0, fill = pal_cns$zebra, inherit.aes = FALSE) +
  annotate("text", x = 0.02, y = 12.02, label = "HR (95% CI)", hjust = 0,
           size = 2.1, fontface = "bold", color = pal_cns$slate, family = cns_font) +
  annotate("text", x = 0.98, y = 12.02, label = "P", hjust = 1,
           size = 2.1, fontface = "bold.italic", color = pal_cns$slate, family = cns_font) +
  geom_text(aes(x = 0.02, label = hr_txt), hjust = 0, size = 2.08,
            color = pal_cns$ink, family = cns_font) +
  geom_text(aes(x = 0.98, label = p_txt, color = p_col,
                fontface = ifelse(sig, "bold", "plain")),
            hjust = 1, size = 2.08, family = cns_font) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim_a, expand = c(0, 0)) +
  theme_void(base_family = cns_font) +
  theme(plot.margin = margin(2, 4, 2, 1))

p6a_body <- (p6a_mat | p6a_for | p6a_tbl) + plot_layout(widths = c(1.22, 1.44, 1.12))
p6a <- wrap_elements(full = p6a_body +
  plot_annotation(
    title = fmt_panel_hdr("A", "Covariate-adjusted <i>PPARG</i> hazard ratios for ER+ overall survival"),
    theme = theme(
      plot.title  = ggtext::element_markdown(hjust = 0, margin = margin(t = 2, b = 2, l = 2)),
      plot.margin = margin(2, 2, 2, 2)
    )
  ))

# ==============================================================================
# PANEL B: 5-Tier Model Discrimination Ladder (C-index +/- SE + delta-AIC/LRT P)
# ==============================================================================
lad_dt <- copy(roc_dt[endpoint == "OS" & cohort %chin% c("TCGA_ERpos", "METABRIC_ERpos")])
lad_dt[, cohort_lbl := factor(coh_strip[cohort], levels = coh_strip)]
lad_dt[, model_lbl  := fcase(
  model == "M1_PPARG_Alone",             "*PPARG* alone",
  model == "M2_Clinical_Base",           "Clinical base",
  model == "M3_Clin_Plus_PPARG",         "Clinical &#43; *PPARG*",
  model == "M4_Clin_Adipose_Base",       "Clin. &#43; adipose",
  model == "M5_Clin_Adipose_Plus_PPARG", "Clin. &#43; adipose &#43; *PPARG*"
)]
mod_lvls <- c("*PPARG* alone", "Clinical base", "Clinical &#43; *PPARG*",
              "Clin. &#43; adipose", "Clin. &#43; adipose &#43; *PPARG*")
lad_dt[, model_lbl := factor(model_lbl, levels = rev(mod_lvls))]
lad_dt[, role := fcase(
  model == "M1_PPARG_Alone", "Single gene",
  grepl("_Base$", model),    "Base",
  default = "Plus"
)]
lad_dt[, pt_col := fcase(
  role == "Single gene", pal_cns$muted,
  role == "Base",        pal_cns$slate,
  cohort == "TCGA_ERpos", pal_cns$pparg,
  default = pal_cns$cxcl8
)]
lad_dt[, pt_fill := ifelse(role == "Plus", pt_col, "white")]

aic_ref <- lad_dt[, .(cohort, model, AIC)]
pair_info <- rbind(
  merge(aic_ref[model == "M3_Clin_Plus_PPARG"], aic_ref[model == "M2_Clinical_Base"],
        by = "cohort", suffixes = c("_plus", "_base"))[, `:=`(comparison = "Clinical",
                                                              model_lbl = "Clinical &#43; *PPARG*")],
  merge(aic_ref[model == "M5_Clin_Adipose_Plus_PPARG"], aic_ref[model == "M4_Clin_Adipose_Base"],
        by = "cohort", suffixes = c("_plus", "_base"))[, `:=`(comparison = "Clinical_Adipose",
                                                              model_lbl = "Clin. &#43; adipose &#43; *PPARG*")]
)
pair_info[, dAIC := AIC_plus - AIC_base]
pair_info <- merge(pair_info, paired[, .(cohort, comparison, LRT_pval)], by = c("cohort", "comparison"))
pair_info[, cohort_lbl := factor(coh_strip[cohort], levels = coh_strip)]
pair_info[, model_lbl  := factor(model_lbl, levels = rev(mod_lvls))]
pair_info[, badge_txt  := fmt_minus(sprintf("\u0394AIC %+.1f \u00b7 P=%.3f", dAIC, LRT_pval))]
pair_info[, badge_col  := ifelse(LRT_pval < 0.05, pal_cns$teal_dark, pal_cns$slate)]
pair_info[, badge_fill := ifelse(LRT_pval < 0.05, pal_cns$teal_soft, "#F1F5F9")]

lad_dt <- merge(lad_dt, pair_info[, .(cohort, model_lbl, badge_txt, badge_col, badge_fill)],
                by = c("cohort", "model_lbl"), all.x = TRUE)
lad_dt[, c_txt := sprintf("%.3f", C_index)]

p6b <- ggplot(lad_dt, aes(x = C_index, y = model_lbl)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = pal_cns$muted, linewidth = 0.42) +
  geom_segment(aes(x = 0.5, xend = C_index, yend = model_lbl, color = pt_col),
               linewidth = 0.65, alpha = 0.35) +
  geom_linerange(aes(xmin = C_index - C_SE, xmax = C_index + C_SE, color = pt_col),
                 linewidth = 0.85) +
  geom_point(aes(color = pt_col, fill = pt_fill), shape = 21, size = 2.4, stroke = 0.82) +
  geom_text(aes(x = C_index + C_SE + 0.008, label = c_txt, color = pt_col,
                fontface = ifelse(role == "Plus", "bold", "plain")),
            hjust = 0, size = 2.08, family = cns_font) +
  geom_label(data = lad_dt[!is.na(badge_txt)],
             aes(x = 0.905, y = model_lbl, label = badge_txt,
                 color = badge_col, fill = badge_fill),
             hjust = 0, size = 1.88, fontface = "bold", linewidth = 0,
             label.r = unit(1.8, "pt"), label.padding = unit(1.8, "pt"),
             family = cns_font, inherit.aes = FALSE) +
  facet_wrap(~ cohort_lbl, ncol = 1) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8), labels = c("0.5", "0.6", "0.7", "0.8")) +
  coord_cartesian(xlim = c(0.49, 0.875), clip = "off") +
  labs(title = fmt_panel_hdr("B", "Model discrimination ladder (C-index \u00b1 SE)"),
       x = "Apparent Harrell C-index (\u00b1 1 SE)", y = NULL) +
  theme_cns(7.8) +
  theme(
    plot.title         = ggtext::element_markdown(hjust = 0, margin = margin(t = 2, b = 5, l = 0)),
    panel.grid.major.x = element_line(color = pal_cns$grid_line, linewidth = 0.30),
    panel.spacing.y    = unit(4.5, "mm"),
    axis.text.y        = ggtext::element_markdown(size = 6.9, color = pal_cns$ink),
    strip.clip         = "off",
    strip.text         = ggtext::element_markdown(face = "bold", size = 7.4, color = pal_cns$ink,
                                                  hjust = 0, margin = margin(t = 2, b = 3)),
    plot.margin        = margin(2, 68, 4, 4)
  )

# ==============================================================================
# PANEL C: Dual-Layer Time-Dependent AUC Trajectories + delta-AUC Diverging Lollipop
# ==============================================================================
roc_top <- rbindlist(lapply(c("Base", "+ PPARG", "PPARG alone"), function(role) {
  if (role == "PPARG alone") {
    sub <- copy(roc_dt[model == "M1_PPARG_Alone" & cohort %chin% c("TCGA_ERpos", "METABRIC_ERpos")])
    rbind(copy(sub)[, comparison_lbl := "Clinical"],
          copy(sub)[, comparison_lbl := "Clinical + adipose"])[, model_role := role]
  } else {
    m <- if (role == "Base") "base_model" else "plus_model"
    keys <- paired[, .(cohort, endpoint, n, events,
                       comparison_lbl = comp_labels[comparison], model = get(m))]
    merge(keys, roc_dt, by = c("cohort", "endpoint", "model", "n", "events"))[, model_role := role]
  }
}), use.names = TRUE, fill = TRUE)

roc_top[, cohort_s := factor(coh_short[cohort], levels = c("TCGA-BRCA", "METABRIC"))]
roc_top[, comparison_lbl := factor(comparison_lbl, levels = c("Clinical", "Clinical + adipose"))]
roc_long <- melt(roc_top, id.vars = c("cohort_s", "comparison_lbl", "model_role"),
                 measure.vars = c("AUC_3yr", "AUC_5yr", "AUC_10yr"),
                 variable.name = "horizon", value.name = "AUC")
roc_long[, years := c(AUC_3yr = 3, AUC_5yr = 5, AUC_10yr = 10)[as.character(horizon)]]
roc_long[, model_role := factor(model_role, levels = c("+ PPARG", "Base", "PPARG alone"))]
roc_long[, grp_id := paste(cohort_s, model_role, sep = "_")]

p6c_top <- ggplot(roc_long, aes(x = years, y = AUC, group = grp_id, color = cohort_s)) +
  geom_line(data = roc_long[model_role == "PPARG alone"],
            linetype = "dotted", linewidth = 0.48, alpha = 0.50) +
  geom_line(data = roc_long[model_role == "Base"],
            linetype = "dashed", linewidth = 0.55, alpha = 0.72) +
  geom_line(data = roc_long[model_role == "+ PPARG"],
            linetype = "solid", linewidth = 0.78) +
  geom_point(data = roc_long[model_role == "Base"],
             shape = 21, fill = "white", size = 1.9, stroke = 0.65) +
  geom_point(data = roc_long[model_role == "+ PPARG"],
             shape = 16, size = 1.95) +
  facet_wrap(~ comparison_lbl, nrow = 1) +
  scale_color_manual(values = coh_pal, name = NULL) +
  scale_x_continuous(breaks = c(3, 5, 10), limits = c(2.0, 11.2)) +
  scale_y_continuous(breaks = c(0.50, 0.60, 0.70, 0.80), limits = c(0.47, 0.845)) +
  labs(x = NULL, y = "Apparent IPCW AUC") +
  theme_cns_grid(7.8) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line.x      = element_line(color = pal_cns$border_soft, linewidth = 0.35),
    strip.text       = ggtext::element_markdown(face = "bold", size = 7.5, color = pal_cns$ink,
                                                hjust = 0.5, margin = margin(t = 1, b = 3)),
    panel.spacing.x  = unit(5.5, "mm"),
    legend.position  = "top",
    legend.margin    = margin(0, 0, -4, 0),
    legend.key.width = unit(5, "mm"),
    plot.margin      = margin(2, 6, 1, 5)
  )

dauc_wide <- dcast(roc_long[model_role %chin% c("Base", "+ PPARG")],
                   cohort_s + comparison_lbl + years ~ model_role, value.var = "AUC")
dauc_wide[, delta_auc := `+ PPARG` - Base]
dauc_wide[, x_shift := ifelse(cohort_s == "TCGA-BRCA", years - 0.38, years + 0.38)]
dauc_wide[, lbl := fmt_minus(sprintf("%+.3f", delta_auc))]
dauc_wide[, lbl_y := ifelse(delta_auc >= 0, delta_auc + 0.0028, delta_auc - 0.0032)]

p6c_bot <- ggplot(dauc_wide, aes(x = x_shift, y = delta_auc, color = cohort_s)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#F5F9FC", alpha = 0.85) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FDF5F4", alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = pal_cns$muted, linewidth = 0.40) +
  geom_segment(aes(xend = x_shift, y = 0, yend = delta_auc), linewidth = 0.80) +
  geom_point(size = 2.05) +
  geom_text(data = dauc_wide[abs(delta_auc) >= 0.005],
            aes(y = lbl_y, label = lbl), size = 1.92, fontface = "bold",
            show.legend = FALSE, family = cns_font) +
  facet_wrap(~ comparison_lbl, nrow = 1) +
  scale_color_manual(values = coh_pal, guide = "none") +
  scale_x_continuous(breaks = c(3, 5, 10), labels = c("3", "5", "10"), limits = c(2.0, 11.2)) +
  scale_y_continuous(breaks = c(-0.01, 0, 0.01, 0.02),
                     labels = c("\u22120.01", "0", "+0.01", "+0.02"),
                     limits = c(-0.013, 0.0265)) +
  labs(x = "OS horizon (years)", y = "\u0394AUC (+ *PPARG*)") +
  theme_cns(7.8) +
  theme(
    strip.text      = element_blank(),
    panel.spacing.x = unit(5.5, "mm"),
    plot.margin     = margin(1, 6, 4, 5)
  )

p6c <- wrap_elements(full = (p6c_top / p6c_bot + plot_layout(heights = c(1.55, 1.0))) +
  plot_annotation(
    title = fmt_panel_hdr("C", "Time-dependent AUC trajectories &amp; incremental &Delta;AUC"),
    theme = theme(
      plot.title  = ggtext::element_markdown(hjust = 0, margin = margin(t = 2, b = 1, l = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

# ==============================================================================
# PANEL D: Free-Scale Bootstrap delta-C Optimism Shrinkage Vector Plot
# ==============================================================================
make_shrink_panel <- function(coh_code, xlim_vec, x_breaks, x_labs,
                              show_xlab = FALSE, show_legend = FALSE) {
  d <- copy(paired[cohort == coh_code])
  d[, comparison_lbl := factor(comp_labels[comparison], levels = rev(comp_labels))]
  col_main <- if (coh_code == "TCGA_ERpos") pal_cns$pparg else pal_cns$cxcl8
  col_band <- if (coh_code == "TCGA_ERpos") pal_cns$pparg_soft else pal_cns$cxcl8_soft
  d[, ci_txt  := fmt_minus(sprintf("%+.5f [%+.5f, %+.5f]", delta_C, delta_C_LCI95, delta_C_UCI95))]
  d[, opt_txt := fmt_minus(sprintf("Opt.-corr: %+.5f", optimism_corrected_delta_C))]
  d[, arr_len := abs(delta_C - optimism_corrected_delta_C)]

  x_text <- xlim_vec[2] + 0.04 * diff(xlim_vec)
  hdr_dt <- data.table(x = x_text, y = 2.48, label = "Apparent \u0394C [95% CI]")

  p <- ggplot(d, aes(y = comparison_lbl)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F8FAFC", alpha = 0.90) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal_cns$muted, linewidth = 0.42) +
    geom_segment(aes(x = delta_C_LCI95, xend = delta_C_UCI95, yend = comparison_lbl),
                 linewidth = 4.2, lineend = "round", color = col_band) +
    geom_segment(aes(x = delta_C_LCI95, xend = delta_C_UCI95, yend = comparison_lbl),
                 linewidth = 0.55, color = col_main, alpha = 0.85) +
    geom_segment(data = d[arr_len > 0.00015],
                 aes(x = delta_C, xend = optimism_corrected_delta_C, yend = comparison_lbl),
                 linewidth = 0.72, color = pal_cns$teal_dark,
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    geom_point(aes(x = delta_C, shape = "Apparent", fill = "Apparent", color = "Apparent"),
               size = 2.35, stroke = 0.65) +
    geom_point(aes(x = optimism_corrected_delta_C, shape = "Optimism-corrected",
                   fill = "Optimism-corrected", color = "Optimism-corrected"),
               size = 2.45, stroke = 0.75) +
    geom_text(aes(x = x_text, y = as.numeric(comparison_lbl) + 0.13, label = ci_txt),
              hjust = 0, size = 2.0, color = pal_cns$ink, family = cns_font) +
    geom_text(aes(x = x_text, y = as.numeric(comparison_lbl) - 0.15, label = opt_txt),
              hjust = 0, size = 1.92, fontface = "bold", color = pal_cns$teal_dark, family = cns_font) +
    scale_shape_manual(values = c("Apparent" = 21, "Optimism-corrected" = 23), name = NULL) +
    scale_color_manual(values = c("Apparent" = col_main, "Optimism-corrected" = pal_cns$teal_dark), name = NULL) +
    scale_fill_manual(values  = c("Apparent" = col_main, "Optimism-corrected" = pal_cns$teal_soft), name = NULL) +
    guides(
      shape = guide_legend(override.aes = list(
        color = c(pal_cns$ink, pal_cns$teal_dark),
        fill  = c(pal_cns$ink, pal_cns$teal_soft)
      )),
      color = "none", fill = "none"
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labs) +
    scale_y_discrete(expand = expansion(add = c(0.55, 0.58))) +
    coord_cartesian(xlim = xlim_vec, clip = "off") +
    labs(subtitle = coh_strip[coh_code],
         x = if (show_xlab) "\u0394C with *PPARG* (95% CI)" else NULL,
         y = NULL) +
    theme_cns(7.8) +
    theme(
      plot.subtitle    = ggtext::element_markdown(face = "bold", size = 7.3, color = pal_cns$ink,
                                                  hjust = 0, margin = margin(t = 1, b = 3, l = 0)),
      axis.text.y      = ggtext::element_markdown(size = 7.1, color = pal_cns$ink),
      legend.position  = if (show_legend) "bottom" else "none",
      legend.margin    = margin(-1, 0, 0, 0),
      legend.key.width = unit(4, "mm"),
      plot.margin      = margin(2, 92, if (show_xlab) 3 else 2, 4)
    )
  if (coh_code == "TCGA_ERpos") {
    p <- p + geom_text(data = hdr_dt, aes(x = x, y = y, label = label),
                       inherit.aes = FALSE, hjust = 0, size = 2.0, fontface = "bold",
                       color = pal_cns$slate, family = cns_font)
  }
  p
}

p6d_tcga <- make_shrink_panel(
  "TCGA_ERpos",
  xlim_vec = c(-0.011, 0.032),
  x_breaks = c(-0.01, 0, 0.01, 0.02, 0.03),
  x_labs   = c("\u22120.01", "0", "0.01", "0.02", "0.03"),
  show_xlab = FALSE, show_legend = FALSE
)
p6d_meta <- make_shrink_panel(
  "METABRIC_ERpos",
  xlim_vec = c(-0.0012, 0.0042),
  x_breaks = c(-0.001, 0, 0.002, 0.004),
  x_labs   = c("\u22120.001", "0", "0.002", "0.004"),
  show_xlab = TRUE, show_legend = TRUE
)

p6d <- wrap_elements(full = (p6d_tcga / p6d_meta + plot_layout(heights = c(1.0, 1.15))) +
  plot_annotation(
    title = fmt_panel_hdr("D", "Bootstrap &Delta;C-index &amp; optimism shrinkage"),
    theme = theme(
      plot.title  = ggtext::element_markdown(hjust = 0, margin = margin(t = 2, b = 1, l = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

# ==============================================================================
# Assemble Final High-Impact 2x2 Composite Figure
# ==============================================================================
row_top <- (p6a | wrap_elements(full = p6b)) + plot_layout(widths = c(1.22, 1.00))
row_bot <- (p6c | p6d) + plot_layout(widths = c(1.04, 1.06))

fig6_cns <- (row_top / row_bot) + plot_layout(heights = c(1.18, 1.04))

save_cns_figure(fig6_cns, "Fig6_ERpos_PPARG_Discrimination", width_mm = 198, height_mm = 214, dpi_png = 350)
cat("High-Impact CNS Fig 6 (Polished) complete.\n")
