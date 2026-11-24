# 02_fig6_pparg_discrim.R
# Figure 6 (was Old Fig 5): ER+ PPARG multivariable Cox and within-cohort incremental discrimination.
# Completely rebuilt to eliminate all title, label, and margin clipping from Old Fig 5.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Final Fig 6 (ER+ PPARG Multivariable Cox & Discrimination) ===\n")

ad_dt  <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
paired <- fread(file.path(res_dir, "planA_erpos_delta_c_bootstrap.csv"))
roc_dt <- fread(file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))

coh_labels  <- c(TCGA_ERpos = "TCGA-BRCA (n=770, 94 deaths)", METABRIC_ERpos = "METABRIC (n=1,127, 626 deaths)")
coh_short   <- c(TCGA_ERpos = "TCGA-BRCA", METABRIC_ERpos = "METABRIC")
comp_labels <- c(Clinical = "Clinical", Clinical_Adipose = "Clinical + adipose")
model_col   <- c("Base" = pal_rev$slate, "+ PPARG" = pal_rev$pparg)

# 6A: ER+ Multivariable Cox HR
ad_sub <- ad_dt[endpoint == "OS" & cohort %chin% names(coh_labels)]
ad_sub[, cohort_lbl := factor(coh_labels[cohort], levels = coh_labels)]
ad_sub[, model_lbl := fcase(
  model == "M1_Univariate", "Univariate PPARG",
  startsWith(model, "M2_"), "+ ADIPOQ & PLIN1",
  model == "M3_Adj_Adipose+Stromal", "+ adipose & stromal",
  model %chin% c("M3_Adj_Age+Stage", "M4_Adj_Age+Stage"), "+ age & stage",
  model == "M4_Adj_Age+Stage+Grade", "+ age, stage & grade",
  startsWith(model, "M5_"), "Fully adjusted"
)]
ad_sub[, model_key := paste(cohort_lbl, model_lbl, sep = "__")]
ad_sub[, model_key := factor(model_key, levels = rev(model_key))]

p6a <- ggplot(ad_sub, aes(x = HR, y = model_key)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.24,
                color = pal_rev$pparg, linewidth = 0.65) +
  geom_point(color = pal_rev$pparg, size = 2.15) +
  facet_wrap(~ cohort_lbl, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.2)) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(title = "A   PPARG and ER+ overall survival", x = "HR per native unit (95% CI)", y = NULL) +
  theme_rev(8.2) +
  theme(strip.text = element_text(size = 7.4))

# 6B: Apparent Harrell C-index (dumbbell with non-overlapping text)
paired[, cohort_lbl := factor(coh_labels[cohort], levels = coh_labels)]
paired[, cohort_s   := factor(coh_short[cohort],  levels = coh_short)]
paired[, comparison_lbl := factor(comp_labels[comparison], levels = rev(comp_labels))]

p6b <- ggplot(paired, aes(y = comparison_lbl)) +
  geom_segment(aes(x = C_index_base, xend = C_index_plus, yend = comparison_lbl),
               linewidth = 0.9, color = pal_rev$border) +
  geom_point(aes(x = C_index_base, color = "Base"), shape = 1, stroke = 0.95, size = 2.8) +
  geom_point(aes(x = C_index_plus, color = "+ PPARG"), size = 2.0) +
  geom_text(aes(x = pmin(C_index_base, C_index_plus) - 0.009,
                label = sprintf("%.4f", C_index_base)), hjust = 1, size = 2.35, color = pal_rev$slate) +
  geom_text(aes(x = pmax(C_index_base, C_index_plus) + 0.009,
                label = sprintf("%.4f", C_index_plus)), hjust = 0, size = 2.35, fontface = "bold", color = pal_rev$pparg) +
  facet_wrap(~ cohort_lbl, ncol = 1) +
  scale_color_manual(values = model_col, breaks = names(model_col), name = NULL) +
  scale_x_continuous(limits = c(0.60, 0.85), breaks = c(0.65, 0.70, 0.75, 0.80)) +
  labs(title = "B   Apparent Harrell C-index", x = "Harrell C-index", y = NULL) +
  theme_rev(8.2) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        strip.text = element_text(size = 7.4))

# 6C: Apparent time-dependent AUC (3, 5, 10 yr)
roc_pairs <- rbindlist(lapply(c("Base", "+ PPARG"), function(which_model) {
  m <- if (which_model == "Base") "base_model" else "plus_model"
  keys <- paired[, .(cohort, endpoint, n, events, cohort_s, comparison_lbl, model = get(m))]
  joined <- merge(keys, roc_dt, by = c("cohort", "endpoint", "model", "n", "events"))
  joined[, model_role := which_model]
  joined
}))
roc_plot <- melt(roc_pairs, id.vars = c("cohort_s", "comparison_lbl", "model_role"),
                 measure.vars = c("AUC_3yr", "AUC_5yr", "AUC_10yr"), variable.name = "horizon", value.name = "AUC")
roc_plot[, years := c(AUC_3yr = 3, AUC_5yr = 5, AUC_10yr = 10)[as.character(horizon)]]
roc_plot[, comparison_lbl := factor(comparison_lbl, levels = comp_labels)]

p6c <- ggplot(roc_plot, aes(x = years, y = AUC, group = model_role, color = model_role)) +
  geom_line(aes(linetype = model_role), linewidth = 0.72) +
  geom_point(aes(shape = model_role), size = 2.1) +
  facet_grid(cohort_s ~ comparison_lbl) +
  scale_x_continuous(breaks = c(3, 5, 10), limits = c(2.2, 10.8)) +
  scale_color_manual(values = model_col, breaks = names(model_col), name = NULL) +
  scale_linetype_manual(values = c("Base" = "dashed", "+ PPARG" = "solid"), breaks = names(model_col), name = NULL) +
  scale_shape_manual(values = c("Base" = 1, "+ PPARG" = 16), breaks = names(model_col), name = NULL) +
  labs(title = "C   Time-dependent AUC", x = "OS horizon (years)", y = "Apparent IPCW AUC") +
  theme_rev_grid(8.2) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        strip.text = element_text(size = 7.3), panel.spacing = unit(2.5, "mm"))

# 6D: Paired bootstrap delta C-index (intervals crossing zero clearly readable)
paired[, dc_txt := sprintf("%+.5f [%+.5f, %+.5f]", delta_C, delta_C_LCI95, delta_C_UCI95)]

p6d <- ggplot(paired, aes(y = comparison_lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
  geom_errorbar(aes(xmin = delta_C_LCI95, xmax = delta_C_UCI95), orientation = "y", width = 0.24,
                linewidth = 0.65, color = pal_rev$pparg) +
  geom_point(aes(x = delta_C, shape = "Apparent"), color = pal_rev$pparg, size = 2.3) +
  geom_point(aes(x = optimism_corrected_delta_C, shape = "Optimism-corrected"),
             color = pal_rev$teal, size = 2.5, stroke = 0.7) +
  geom_text(aes(x = 0.033, label = dc_txt), hjust = 0, size = 2.15, color = pal_rev$ink) +
  facet_wrap(~ cohort_lbl, ncol = 1) +
  scale_shape_manual(values = c("Apparent" = 16, "Optimism-corrected" = 23), name = NULL) +
  scale_x_continuous(limits = c(-0.014, 0.088), breaks = c(-0.01, 0.01, 0.03),
                     labels = c("-0.01", "0.01", "0.03")) +
  labs(title = "D   Paired bootstrap \u0394C-index", x = "\u0394C with PPARG (95% CI)", y = NULL) +
  theme_rev(8.2) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        strip.text = element_text(size = 7.4))

fig6_all <- (wrap_elements(full = p6a) | wrap_elements(full = p6b)) /
            (wrap_elements(full = p6c) | wrap_elements(full = p6d)) +
  plot_layout(heights = c(1.15, 1.05), widths = c(0.95, 1.15))

save_rev_figure(fig6_all, "Fig6_ERpos_PPARG_Discrimination", width_mm = 180, height_mm = 184, dpi_png = 350)
cat("Fig 6 complete.\n")
