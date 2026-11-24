# 25c_generate_cns_figs4to5.R
# CNS 级顶刊版 Fig 4 (Raincloud + 显著性括号 + 95% CI 阴影带 KM + Meta) & Fig 5 (矢量渐变列线图 + 哑铃增量图)

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/25_0_cns_plot_helpers.R"))

# ==============================================================================
# FIG 4: 临床表达分化 Raincloud (4A/4B) + 95% CI KM (4C–4F) + 跨队列 Meta (4G)
# ==============================================================================
cat("=== Building Fig 4 (Raincloud + 95% CI KM + Pooled Meta) ===\n")

tc  <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb  <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20 <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42 <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
diff_dt <- fread(file.path(res_dir, "module4_differential_expression.csv"))
s13     <- fread(file.path(res_dir, "planA_pam50_subtypes_stats.csv"))

# 4A: Tumor vs Normal Raincloud
df_tn <- rbind(
  data.table(Cohort = "GSE42568",
             G_simple = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal", "Tumor"),
             Expr = as.numeric(g42$expr_target["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Normal", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Tumor",  Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)
df_tn[, G_simple := factor(G_simple, levels = c("Normal", "Tumor"))]
vio_4a <- build_half_violin_dt(df_tn, "G_simple", "Expr", facet_col = "Cohort", width = 0.34, side = "right")

p4a_gse_p <- diff_dt[cohort == "GSE42568"  & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
p4a_tcg_p <- diff_dt[cohort == "TCGA_BRCA" & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
brk_4a <- df_tn[, .(ymax = max(Expr, na.rm = TRUE), ymin = min(Expr, na.rm = TRUE)), by = Cohort]
brk_4a[, `:=`(y_brk = ymax + 0.08 * (ymax - ymin), y_tick = ymax + 0.03 * (ymax - ymin),
              y_txt = ymax + 0.14 * (ymax - ymin),
              p_expr = c(format_pval_expr(p4a_gse_p, "Wilcoxon~italic(P)"),
                         format_pval_expr(p4a_tcg_p, "Wilcoxon~italic(P)")))]

p4a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple, color = G_simple)) +
  geom_polygon(data = vio_4a, aes(x = vx, y = vy, group = grp_id, fill = G_simple),
               inherit.aes = FALSE, alpha = 0.35, color = "#475569", linewidth = 0.35) +
  ggrastr::rasterise(geom_jitter(aes(x = as.numeric(G_simple) - 0.18), width = 0.07, size = 0.55, alpha = 0.32), dpi = 300) +
  geom_boxplot(aes(x = as.numeric(G_simple) - 0.18), width = 0.14, outlier.shape = NA, alpha = 0.88, color = "#0F172A", linewidth = 0.45) +
  geom_segment(data = brk_4a, aes(x = 1, xend = 2, y = y_brk, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.5) +
  geom_segment(data = brk_4a, aes(x = 1, xend = 1, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.5) +
  geom_segment(data = brk_4a, aes(x = 2, xend = 2, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.5) +
  geom_text(data = brk_4a, aes(x = 1.5, y = y_txt, label = p_expr), inherit.aes = FALSE, parse = TRUE, size = 2.65, color = "#0F172A") +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_x_discrete(limits = c("Normal", "Tumor")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  scale_fill_manual(values = c("Normal" = "#38BDF8", "Tumor" = pal_cns[["red"]])) +
  scale_color_manual(values = c("Normal" = "#0284C7", "Tumor" = pal_cns[["red"]])) +
  labs(title = "A. PPARG Expression (Tumor vs. Normal)", x = NULL, y = expression(bold(Normalized~log[2]~Expression))) +
  theme_cns(9.3) + theme(legend.position = "none")

# 4B: PAM50 Subtypes Raincloud + KW & Basal vs LumA P
pam_tc <- rbind(
  data.table(Gene = "CXCL8", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ])),
  data.table(Gene = "PPARG", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)[Subtype %in% c("Normal","LumA","LumB","Her2","Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal","LumA","LumB","Her2","Basal"),
                           labels = c("Normal","LumA","LumB","HER2+","Basal"))]
pam_tc[, Gene := factor(Gene, levels = c("CXCL8", "PPARG"))]
vio_4b <- build_half_violin_dt(pam_tc, "Subtype", "Expr", facet_col = "Gene", width = 0.36, side = "right")

kw_cx <- s13[cohort == "TCGA_BRCA" & gene_symbol == "CXCL8"]$kw_pval[1]
w_cx  <- wilcox.test(pam_tc[Gene == "CXCL8" & Subtype == "Basal"]$Expr,
                     pam_tc[Gene == "CXCL8" & Subtype == "LumA"]$Expr)$p.value
kw_pp <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$kw_pval[1]
w_pp  <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$wilcox_Basal_vs_LumA[1]

ann_4b <- data.table(
  Gene = factor(c("CXCL8", "PPARG"), levels = c("CXCL8", "PPARG")),
  ymax = c(max(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), max(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  ymin = c(min(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), min(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  kw_lbl = c(format_pval_expr(kw_cx, "Kruskal-Wallis~italic(P)"), format_pval_expr(kw_pp, "Kruskal-Wallis~italic(P)")),
  w_lbl  = c(format_pval_expr(w_cx, "Basal~vs.~LumA~italic(P)"),   format_pval_expr(w_pp, "Basal~vs.~LumA~italic(P)"))
)
ann_4b[, `:=`(y_brk = ymax + 0.08 * (ymax - ymin), y_tick = ymax + 0.03 * (ymax - ymin),
              y_w = ymax + 0.16 * (ymax - ymin), y_kw = ymax + 0.28 * (ymax - ymin))]

p4b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_polygon(data = vio_4b, aes(x = vx, y = vy, group = grp_id, fill = Subtype),
               inherit.aes = FALSE, alpha = 0.35, color = "#475569", linewidth = 0.30) +
  geom_boxplot(aes(x = as.numeric(Subtype) - 0.14), width = 0.22, outlier.size = 0.4, alpha = 0.88, color = "#0F172A", linewidth = 0.42) +
  geom_segment(data = ann_4b, aes(x = 2, xend = 5, y = y_brk, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.48) +
  geom_segment(data = ann_4b, aes(x = 2, xend = 2, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.48) +
  geom_segment(data = ann_4b, aes(x = 5, xend = 5, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.48) +
  geom_text(data = ann_4b, aes(x = 3.5, y = y_w, label = w_lbl), inherit.aes = FALSE, parse = TRUE, size = 2.45, color = "#0F172A") +
  geom_text(data = ann_4b, aes(x = 3.0, y = y_kw, label = kw_lbl), inherit.aes = FALSE, parse = TRUE, size = 2.55, color = "#334155") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_x_discrete(limits = c("Normal","LumA","LumB","HER2+","Basal")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  scale_fill_manual(values = c("Normal"="#94A3B8","LumA"=pal_cns[["blue"]],"LumB"="#6366F1","HER2+"=pal_cns[["amber"]],"Basal"=pal_cns[["red"]])) +
  labs(title = "B. PAM50 Subtype Expression (TCGA-BRCA)", x = NULL, y = expression(bold(log[2](RSEM + 1)))) +
  theme_cns(9.3) + theme(legend.position = "none")

# 4C-4F: 95% CI Ribbon + HR + Log-rank P + Balanced Risk Table KM Curves
idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp],
                      expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4c <- plot_km_cns(
  df_tc_p, "time", "status", "Group",
  "C. TCGA-BRCA ER+ OS (PPARG)", "Months", "Overall Survival", "PPARG",
  c("Low", "High"), c(pal_cns[["red"]], pal_cns[["blue"]]),
  pval_coord = c(8, 0.14), break_by = 50, max_x = 250
)

idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp],
                      expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4d <- plot_km_cns(
  df_mb_p, "time", "status", "Group",
  "D. METABRIC ER+ OS (PPARG)", "Months", "Overall Survival", "PPARG",
  c("Low", "High"), c(pal_cns[["red"]], pal_cns[["blue"]]),
  pval_coord = c(10, 0.14), break_by = 100, max_x = 300
)

pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4e <- plot_km_cns(
  df_g20_c, "time", "status", "Group",
  "E. GSE20685 OS (CXCL8)", "Years", "Overall Survival", "CXCL8",
  c("Low", "High"), c(pal_cns[["blue"]], pal_cns[["red"]]),
  pval_coord = c(0.5, 0.16), break_by = 3, max_x = 15
)

idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern],
                      expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4f <- plot_km_cns(
  df_mb_c, "time", "status", "Group",
  "F. METABRIC ER\u2212 OS (CXCL8)", "Months", "Overall Survival", "CXCL8",
  c("Low", "High"), c(pal_cns[["blue"]], pal_cns[["red"]]),
  pval_coord = c(10, 0.14), break_by = 100, max_x = 300
)
rm(tc, mb, g20, g42, cl_tc); gc()

# Figure 4G: Show the standardized OS model. Endpoint-specific results are tabulated separately.
meta_hr <- fread(file.path(res_dir, "planA_standardized_os_meta.csv"))
hr_sub <- meta_hr[(gene_symbol == "CXCL8" & stratum %chin% c("Overall", "ER_neg")) |
                  (gene_symbol == "PPARG" & stratum %chin% c("Overall", "ER_pos"))]
hr_sub[, label := sprintf("%s / %s (%d cohorts)", gene_symbol, clean_subtype_str(stratum), k)]
hr_sub[, label := factor(label, levels = rev(label))]
p4g <- ggplot(hr_sub, aes(x = RE_HR, y = label, color = gene_symbol)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8") +
  geom_errorbar(aes(xmin = RE_LCI, xmax = RE_UCI), orientation = "y", width = 0.25, linewidth = 0.65) +
  geom_point(size = 2.3) +
  scale_color_manual(values = c(CXCL8 = pal_cns[["red_dark"]], PPARG = pal_cns[["blue_dark"]])) +
  labs(title = "G. Random-effects meta-analysis of overall survival",
       x = "Hazard ratio per SD of expression (95% Wald CI)", y = NULL) +
  theme_cns(10) + theme(legend.position = "none")

fig4_r1 <- patchwork::wrap_elements(full = (p4a | p4b))
fig4_all <- fig4_r1 / (km_4c | km_4d) / (km_4e | km_4f) / p4g + plot_layout(heights = c(1, 1.16, 1.16, 0.75))
save_cns(fig4_all, "Fig4_Clinical_Expression_and_Prognosis", w = 180 / 25.4, h = 252 / 25.4, dpi_png = 450)
rm(fig4_all, fig4_r1, p4a, p4b, km_4c, km_4d, km_4e, km_4f, p4g); gc()

# ==============================================================================
# FIG 5: ER+ PPARG associations and within-cohort discrimination
# ==============================================================================
cat("=== Building Fig 5 (ER+ OS associations and discrimination) ===\n")

ad_dt <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
paired <- fread(file.path(res_dir, "planA_erpos_delta_c_bootstrap.csv"))
roc_dt <- fread(file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))
ph <- fread(file.path(res_dir, "planA_erpos_os_ph_diagnostics.csv"))
coh_labels <- c(TCGA_ERpos = "TCGA-BRCA", METABRIC_ERpos = "METABRIC")
comp_labels <- c(Clinical = "Clinical", Clinical_Adipose = "Clinical + adipose covariates")
model_col <- c("Base" = "#64748B", "+ PPARG" = pal_cns[["blue"]])

ad_sub <- ad_dt[endpoint == "OS" & cohort %chin% names(coh_labels)]
ad_sub[, cohort_lbl := factor(coh_labels[cohort], levels = coh_labels)]
ad_sub[, model_lbl := fcase(
  model == "M1_Univariate", "Univariate PPARG",
  startsWith(model, "M2_"), "+ ADIPOQ and PLIN1",
  model == "M3_Adj_Adipose+Stromal", "+ adipose and stromal score",
  model %chin% c("M3_Adj_Age+Stage", "M4_Adj_Age+Stage"), "+ age and stage",
  model == "M4_Adj_Age+Stage+Grade", "+ age, stage and grade",
  startsWith(model, "M5_"), "Fully adjusted")]
ad_sub[, model_key := paste(cohort_lbl, model_lbl, sep = "__")]
ad_sub[, model_key := factor(model_key, levels = rev(model_key))]
stopifnot(!anyNA(ad_sub$model_lbl), all(ad_sub$sample_set == "ERpos_OS_common_complete_cases"))
p5a <- ggplot(ad_sub, aes(x = HR, y = model_key)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.24,
                color = pal_cns[["blue"]], linewidth = 0.75) +
  geom_point(color = pal_cns[["blue"]], size = 2.55) +
  facet_wrap(~ cohort_lbl, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.2)) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(title = "A. PPARG and overall survival in ER+ disease", x = "Hazard ratio per native expression unit (95% CI)",
       y = NULL) + theme_cns(10.5)

paired[, cohort_lbl := factor(coh_labels[cohort], levels = coh_labels)]
paired[, comparison_lbl := factor(comp_labels[comparison], levels = rev(comp_labels))]
paired[, cohort_count := sprintf("%s  |  n = %s, OS events = %s", cohort_lbl, n, events)]
paired[, cohort_count := factor(cohort_count, levels = unique(cohort_count[order(cohort_lbl)]))]
stopifnot(!anyNA(paired$comparison_lbl), !anyNA(paired$cohort_lbl),
          !anyDuplicated(paired[, .(cohort, comparison)]),
          isTRUE(all.equal(paired$C_index_plus - paired$C_index_base, paired$delta_C, tolerance = 1e-10)),
          nrow(merge(unique(ad_sub[, .(cohort, n, events)]), unique(paired[, .(cohort, n, events)]))) == length(coh_labels))

# B 和 D 共用同一张配对自助法汇总表及相同的病例数。
ranges <- diff(range(c(paired$C_index_base, paired$C_index_plus)))
p5b <- ggplot(paired, aes(y = comparison_lbl)) +
  geom_segment(aes(x = C_index_base, xend = C_index_plus, yend = comparison_lbl),
               linewidth = 1.0, color = "#CBD5E1") +
  geom_point(aes(x = C_index_base, color = "Base"), shape = 1, stroke = 1.0, size = 3.4) +
  geom_point(aes(x = C_index_plus, color = "+ PPARG"), size = 2.25) +
  geom_text(aes(x = pmin(C_index_base, C_index_plus) - 0.035 * ranges,
                label = sprintf("%.4f", C_index_base)), hjust = 1, size = 3.1, color = model_col[["Base"]]) +
  geom_text(aes(x = pmax(C_index_base, C_index_plus) + 0.035 * ranges,
                label = sprintf("%.4f", C_index_plus)), hjust = 0, size = 3.1, color = model_col[["+ PPARG"]]) +
  facet_wrap(~ cohort_count, ncol = 1) + scale_color_manual(values = model_col, breaks = names(model_col)) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "B. Apparent Harrell C-index", x = "Harrell C-index for overall survival", y = NULL, color = NULL) +
  theme_cns(10.5) + theme(legend.position = "bottom")

roc_pairs <- rbindlist(lapply(c("Base", "+ PPARG"), function(which_model) {
  m <- if (which_model == "Base") "base_model" else "plus_model"
  keys <- paired[, .(cohort, endpoint, n, events, cohort_lbl, comparison_lbl, model = get(m))]
  joined <- merge(keys, roc_dt, by = c("cohort", "endpoint", "model", "n", "events"))
  stopifnot(nrow(joined) == nrow(keys))
  joined[, model_role := which_model]
  joined
}))
roc_plot <- melt(roc_pairs, id.vars = c("cohort_lbl", "comparison_lbl", "model_role"),
                 measure.vars = c("AUC_3yr", "AUC_5yr", "AUC_10yr"), variable.name = "horizon", value.name = "AUC")
stopifnot(all(is.finite(roc_plot$AUC)))
roc_plot[, years := c(AUC_3yr = 3, AUC_5yr = 5, AUC_10yr = 10)[as.character(horizon)]]
roc_plot[, comparison_lbl := factor(comparison_lbl, levels = comp_labels)]
p5c <- ggplot(roc_plot, aes(x = years, y = AUC, group = model_role, color = model_role)) +
  geom_line(aes(linetype = model_role), linewidth = 0.7) +
  geom_point(aes(shape = model_role), size = 2.45) +
  facet_grid(cohort_lbl ~ comparison_lbl) +
  scale_x_continuous(breaks = c(3, 5, 10)) +
  scale_color_manual(values = model_col, breaks = names(model_col)) +
  scale_linetype_manual(values = c("Base" = "dashed", "+ PPARG" = "solid")) +
  scale_shape_manual(values = c("Base" = 1, "+ PPARG" = 16)) +
  labs(title = "C. Apparent time-dependent AUC", x = "Overall survival horizon (years)", y = "IPCW AUC") +
  theme_cns(10.5) + theme(legend.position = "none")

range_d <- range(c(paired$delta_C_LCI95, paired$delta_C_UCI95, paired$optimism_corrected_delta_C, 0))
caption_d <- sprintf("Paired patient bootstrap (B = %s). Intervals describe the apparent change.\nMETABRIC clinical models showed global PH violations.",
                     format(unique(paired$bootstrap_reps), big.mark = ","))
stopifnot(uniqueN(paired$bootstrap_reps) == 1L,
          nrow(ph[cohort == "METABRIC_ERpos" & term == "GLOBAL"]) == uniqueN(roc_dt$model))
p5d <- ggplot(paired, aes(y = comparison_lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = delta_C_LCI95, xmax = delta_C_UCI95), orientation = "y", width = 0.24,
                linewidth = 0.75, color = pal_cns[["blue"]]) +
  geom_point(aes(x = delta_C, shape = "Apparent"), color = pal_cns[["blue"]], size = 2.75) +
  geom_point(aes(x = optimism_corrected_delta_C, shape = "Optimism-corrected"),
             color = pal_cns[["green"]], size = 3.1) +
  geom_text(aes(x = range_d[2] + 0.12 * diff(range_d),
                label = sprintf("%+.5f [%+.5f, %+.5f]", delta_C, delta_C_LCI95, delta_C_UCI95)),
            hjust = 0, size = 3.0, color = "#334155") +
  facet_wrap(~ cohort_count, ncol = 1) +
  scale_shape_manual(values = c("Apparent" = 16, "Optimism-corrected" = 23)) +
  scale_x_continuous(expand = expansion(mult = c(0.07, 0.55))) +
  labs(title = "D. Paired change in Harrell C-index", x = "Change with PPARG (95% CI for apparent change)",
       y = NULL, shape = NULL, caption = caption_d) +
  theme_cns(10.5) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9))

fig5_all <- (p5a | p5b) / (p5c | p5d) + plot_layout(heights = c(1.2, 1))
save_cns(fig5_all, "Fig5_ERpos_Prognostic_Association_and_Discrimination", w = 180 / 25.4, h = 192 / 25.4, dpi_png = 450)
cat("Fig 5 complete.\n")
