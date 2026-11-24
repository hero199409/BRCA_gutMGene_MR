# 26b_generate_cns_figs8toS3.R
# Frontiers 图版，使用本地结果表及绘图缓存，原输出文件名保持一致。
plot_script_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts")
source(file.path(plot_script_dir, "25_0_cns_plot_helpers.R"))

save_frontiers_26b <- function(p, name, height_mm) {
  save_cns(p, name, w = 180 / 25.4, h = height_mm / 25.4, dpi_png = 450)
}

# 文本与置信区间分列，避免标注占用或截断效应量坐标轴。
forest_26b <- function(d, title, x_title, null, stats_header, widths,
                       log_x = FALSE, x_range = NULL, subtitle = NULL, compact = FALSE) {
  d <- copy(d)
  d[, row := rev(seq_len(.N))]
  zebra <- d[row %% 2L == 1L]
  y_scale <- scale_y_continuous(limits = c(0.45, nrow(d) + 0.55), expand = c(0, 0))
  row_bands <- geom_rect(data = zebra, aes(ymin = row - 0.49, ymax = row + 0.49),
                         xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#F1F5F9")
  table_side <- function(column, header) {
    ggplot(d, aes(y = row)) + row_bands +
      geom_text(aes(x = 0.02, label = .data[[column]]), hjust = 0, size = 2.80,
                lineheight = 0.95, color = "#0F172A") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + y_scale +
      labs(title = header) + theme_void(base_family = "sans") +
      theme(plot.title = element_text(size = 8.0, face = "bold", margin = margin(b = 3)),
            plot.margin = margin(1, 2, 1, 2))
  }
  center <- ggplot(d, aes(x = effect, y = row)) + row_bands +
    geom_vline(xintercept = null, color = "#64748B", linetype = "dashed", linewidth = 0.45) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.24,
                  linewidth = 0.55, color = "#334155") +
    geom_point(size = 1.75, color = pal_cns[["blue"]]) + y_scale +
    labs(title = " ", x = x_title, y = NULL) + theme_cns(8.8) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 8.0),
          plot.title = element_text(size = 8.0, margin = margin(b = 3)),
          plot.margin = margin(1, 3, 1, 3))
  if (log_x) {
    center <- center + scale_x_log10(breaks = if (compact) c(0.5, 0.75, 1, 1.5, 2, 3) else scales::breaks_log(n = 6))
  }
  if (!is.null(x_range)) center <- center + coord_cartesian(xlim = x_range)
  result <- (table_side("row_label", " ") | center | table_side("stats", stats_header)) +
    plot_layout(widths = widths) +
    plot_annotation(title = title, subtitle = subtitle, theme = theme(
      plot.title = element_text(size = 10, face = "bold", color = "#0F172A", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 8.0, color = "#334155", margin = margin(b = 3)),
      plot.margin = margin(3, 4, 3, 4)))
  patchwork::wrap_elements(full = result)
}

# ==============================================================================
# FIG 8: 真实临床治疗响应 (8A 表图合一 + 8B/8C + 8D 蜂群Raincloud + 8E/8F)
# ==============================================================================
cat("=== Building Fig 8 (Table-Aligned Forest + Beeswarm Raincloud + GDSC2) ===\n")

th_dt <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))
th_pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

# 8A：连续表达 Cox 效应量。分组 KM 统计量留在对应的 KM 面板。
endo_dt <- copy(th_dt[gene == "PPARG" & cohort %in% c("METABRIC", "GSE9195", "GSE17705") & stratum == "ER_pos"])
endo_dt[, arm_clean := fcase(treatment_arm == "Hormone_Therapy_YES", "Endocrine (+)",
                            treatment_arm == "Hormone_Therapy_NO", "Endocrine (\u2212)",
                            treatment_arm == "Tamoxifen_5yr", "Tamoxifen", default = NA_character_)]
endo_dt[is.na(arm_clean), arm_clean := treatment_arm]
endo_dt[, `:=`(row_label = sprintf("%s | %s | %s | n=%s", cohort, arm_clean, endpoint,
                                   format(n_samples, big.mark = ",", trim = TRUE)),
               effect = estimate, low = lci95, high = uci95,
               stats = sprintf("%.2f [%.2f, %.2f] | %.3g", estimate, lci95, uci95, pval))]
p8a <- forest_26b(endo_dt, "A. PPARG and endocrine therapy outcomes (ER+)", "HR (log scale)", 1,
                  "HR [95% CI] | Cox P", c(2.05, 1.45, 1.08), log_x = TRUE)

# 8B：单变量连续表达的 OR 和 AUC，P 值来自同一 logistic 模型。
nact_pcr_dt <- copy(th_dt[gene == "CXCL8" & cohort %in% c("GSE25066", "GSE32646") &
                          endpoint == "pCR_vs_RD" & stratum %in% c("Overall", "ER_neg", "TNBC")])
nact_pcr_dt[, `:=`(row_label = sprintf("%s %s\nn=%d (%d pCR)", cohort, clean_subtype_str(stratum), n_samples, n_events),
                   effect = estimate, low = lci95, high = uci95,
                   stats = sprintf("%.2f | %.3g", auc, pval))]
p8b <- forest_26b(nact_pcr_dt, "B. CXCL8 and pCR after NACT", "OR [95% CI] (log scale)", 1,
                  "AUC | P", c(1.10, 1.28, 0.66), log_x = TRUE, compact = TRUE)

# 8C: GSE25066 NACT DRFS 95% CI KM
g25_km <- copy(th_pc$gse25066)[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs)]
g25_km[, Group := ifelse(CXCL8 >= median(CXCL8), "High", "Low")]
km_8c <- plot_km_cns(
  g25_km, "t_drfs", "e_drfs", "Group",
  "C. GSE25066 NACT DRFS (CXCL8)", "Years", "Distant Relapse-Free Survival", "CXCL8",
  c("Low", "High"), c(pal_cns[["blue"]], pal_cns[["red"]]),
  pval_coord = c(0.3, 0.11), break_by = 2, max_x = 8, base_size = 8.6
)

# 8D1 / 8D2：图样本、臂内 logistic P 与交互检验按同一亚组及基因匹配。
th_int <- fread(file.path(res_dir, "planA_therapy_interactions_and_multivar.csv"))
arm_map <- data.table(
  Arm = c("control", "durvalumab/olaparib"),
  treatment_arm = c("Control_Paclitaxel_Arm", "Durvalumab_Olaparib_Arm"),
  label = c("Control", "Durva + Ola")
)
make_ispy_panel <- function(gene_name, population, interaction_suffix, title, model_label, palette) {
  d <- copy(as.data.table(th_pc$gse173839))[pCR %in% c(0L, 1L) & is.finite(get(gene_name))]
  if (population == "TNBC_HRneg") d <- d[HR_neg == TRUE]
  tests <- th_dt[cohort == "GSE173839" & domain == "Neoadjuvant_pCR" &
                  stratum == population & gene == gene_name & treatment_arm %in% arm_map$treatment_arm]
  if (nrow(tests) != nrow(arm_map) || anyDuplicated(tests$treatment_arm)) {
    stop("Missing or ambiguous I-SPY2 within-arm results for ", gene_name, " / ", population)
  }
  tests <- tests[match(arm_map$treatment_arm, treatment_arm)]
  int_test <- paste0("pCR_Interaction_", gene_name, "_x_DurvaArm(", interaction_suffix, ")")
  interaction <- th_int[cohort == "GSE173839_ISPY2" & stratum == population & test_name == int_test]
  counts <- d[, .(n = .N, events = sum(pCR)), by = Arm][match(arm_map$Arm, Arm)]
  if (nrow(interaction) != 1L || nrow(tests) != 2L || anyNA(counts) ||
      !all(is.finite(c(tests$pval, interaction$pval, interaction$fdr_bh, interaction$fdr_family_n))) ||
      !all(counts$n == tests$n_samples & counts$events == tests$n_events) ||
      !identical(as.integer(counts$n), as.integer(c(interaction$n_control, interaction$n_durva))) ||
      !identical(as.integer(counts$events), as.integer(c(interaction$events_control, interaction$events_durva)))) {
    stop("I-SPY2 plotting cache and model summaries do not match for ", gene_name, " / ", population)
  }
  facet_labels <- sprintf("%s (n=%d)", arm_map$label, counts$n)
  d[, `:=`(Arm_lbl = factor(Arm, levels = arm_map$Arm, labels = facet_labels),
           Resp = factor(pCR, levels = c(0, 1), labels = c("RD", "pCR")), expr = get(gene_name))]
  violin <- build_half_violin_dt(d, "Resp", "expr", facet_col = "Arm_lbl", width = 0.27, side = "right")
  y_range <- range(c(d$expr, violin$vy))
  y_span <- diff(y_range)
  bracket <- data.table(Arm_lbl = factor(facet_labels, levels = facet_labels),
                        y = y_range[2] + 0.05 * y_span, tip = y_range[2] + 0.02 * y_span,
                        y_text = y_range[2] + 0.12 * y_span,
                        label = sprintf("Arm P = %.3g", tests$pval))
  ggplot(d, aes(y = expr, fill = Resp, color = Resp)) +
    geom_polygon(data = violin, aes(x = vx, y = vy, group = grp_id, fill = Resp),
                 inherit.aes = FALSE, alpha = 0.32, color = "#475569", linewidth = 0.3) +
    ggbeeswarm::geom_quasirandom(aes(x = as.numeric(Resp) - 0.12, group = Resp),
                                width = 0.08, size = 1.08, alpha = 0.67) +
    geom_boxplot(aes(x = as.numeric(Resp) - 0.12, group = Resp), width = 0.16,
                 outlier.shape = NA, alpha = 0.62, color = "#0F172A", linewidth = 0.42) +
    geom_segment(data = bracket, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket, aes(x = 1, xend = 1, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket, aes(x = 2, xend = 2, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = bracket, aes(x = 1.5, y = y_text, label = label), inherit.aes = FALSE, size = 2.8) +
    facet_wrap(~ Arm_lbl) +
    scale_x_continuous(breaks = c(1, 2), labels = c("RD", "pCR"), limits = c(0.67, 2.36)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    scale_fill_manual(values = c("RD" = pal_cns[["null"]], "pCR" = pal_cns[[palette]])) +
    scale_color_manual(values = c("RD" = "#475569", "pCR" = pal_cns[[paste0(palette, "_dark")]])) +
    labs(title = title, x = NULL, y = paste(gene_name, "expression"),
         subtitle = sprintf("%s\nP = %.3g | BH q = %.3g (m=%d)", model_label,
                            interaction$pval, interaction$fdr_bh, interaction$fdr_family_n)) +
    theme_cns(8.8) +
    theme(plot.subtitle = element_text(size = 8.0, color = "#334155", margin = margin(b = 4)),
          strip.text = element_text(size = 8.0, face = "bold"),
          panel.spacing = unit(3, "mm"), legend.position = "none")
}

p8d1 <- make_ispy_panel("CXCL8", "HER2neg_All", "Adj_HR",
                       "D1. I-SPY2 CXCL8 (HER2\u2212)", "Treatment interaction (HR-adjusted)", "red")
p8d2 <- make_ispy_panel("PPARG", "TNBC_HRneg", "HRneg_Subset",
                       "D2. I-SPY2 PPARG (HER2\u2212 / HR\u2212)", "Treatment interaction (exploratory subgroup)", "blue")

# 8E：按药物 ID 保留全家族。用调整模型 P 排序，BH 在筛选展示条目前核对。
sub_d <- fread(file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
adj <- sub_d[gene_symbol %in% c("CXCL8", "PPARG") & is.finite(subtype_adjusted_pval)]
if (anyDuplicated(adj, by = c("gene_symbol", "DRUG_ID"))) stop("Non-unique adjusted GDSC2 drug IDs")
bh_check <- adj[, .(matches = isTRUE(all.equal(subtype_adjusted_FDR,
                                               p.adjust(subtype_adjusted_pval, method = "BH"), tolerance = 1e-10)) &&
                               all(subtype_adjusted_family_n == .N)), by = gene_symbol]
if (nrow(bh_check) != 2L || !isTRUE(all(bh_check$matches))) stop("GDSC2 adjusted BH family does not match the full table")

make_drug_forest <- function(gene, example_id, x_range) {
  gene_data <- adj[gene_symbol == gene]
  key <- gene_data[is_key_drug == TRUE & DRUG_ID != example_id][order(subtype_adjusted_pval, DRUG_ID)]
  example <- gene_data[DRUG_ID == example_id]
  if (nrow(example) != 1L || nrow(key) < 5L) stop("Missing adjusted GDSC2 drug estimates for ", gene)
  d <- rbind(key[seq_len(5L)], example)
  setorder(d, subtype_adjusted_pval, DRUG_ID)
  d[, row_label := gsub("_", "-", DRUG_NAME, fixed = TRUE)]
  d[duplicated(row_label) | duplicated(row_label, fromLast = TRUE),
    row_label := sprintf("%s (%s)", row_label, DRUG_ID)]
  d[, `:=`(effect = subtype_adjusted_beta, low = subtype_adjusted_ci_low, high = subtype_adjusted_ci_high,
           stats = sprintf("%.3f", subtype_adjusted_FDR))]
  if (!all(is.finite(c(d$effect, d$low, d$high)))) stop("Missing adjusted GDSC2 confidence interval")
  forest_26b(d, sprintf("%s (m=%d)", gene, nrow(gene_data)), "Adjusted \u03b2 [95% CI]",
             0, "BH q", c(1.09, 1.45, 0.40), x_range = x_range)
}
adj_range <- adj[is_key_drug == TRUE, range(c(subtype_adjusted_ci_low, subtype_adjusted_ci_high), finite = TRUE)]
p8e <- patchwork::wrap_elements(full = (make_drug_forest("CXCL8", 1079L, adj_range) |
                                       make_drug_forest("PPARG", 1200L, adj_range)) +
  plot_annotation(title = "E. GDSC2 drug response adjusted for subtype",
                  subtitle = "The F example and the five other key drugs with lowest adjusted P for each gene.\nBH uses all adjusted drug tests per gene. Outcome: ln(IC50).",
                  theme = theme(plot.title = element_text(size = 10, face = "bold", margin = margin(b = 3)),
                                plot.subtitle = element_text(size = 8.0, color = "#334155"),
                                plot.margin = margin(3, 4, 0, 4))))

gb_mat <- readRDS(file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))
gb_mat[, Subtype := factor(Subtype_Group, levels = c("ER_pos", "HER2_pos", "TNBC"), labels = c("ER+", "HER2+", "TNBC"))]

# F 面板使用汇总表中的未调整 Spearman，散点取对应药物 ID。
spearman_label_26b <- function(gene_name, drug_id, d, target_label) {
  rho <- sub_d[gene_symbol == gene_name & DRUG_ID == drug_id]
  d <- d[is.finite(get(gene_name)) & is.finite(LN_IC50)]
  if (nrow(rho) != 1L || nrow(d) != rho$n_all ||
      nrow(d[Subtype_Group == rho$target_subtype]) != rho$n_target_subtype) {
    stop("GDSC2 scatter cache and summary do not match for ", gene_name, " / ", drug_id)
  }
  sprintf("Spearman (unadjusted)\nAll lines: \u03c1 = %.3f (P = %.3f)\n%s: \u03c1 = %.3f (P = %.3f)",
          rho$all_spearman_rho, rho$all_spearman_pval,
          target_label, rho$within_subtype_rho, rho$within_subtype_pval)
}
d_f1 <- gb_mat[DRUG_ID == 1079L]
f1_lbl <- spearman_label_26b("CXCL8", 1079L, d_f1, "TNBC")

p8f1 <- ggplot(d_f1, aes(x = CXCL8, y = LN_IC50, color = Subtype)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75, alpha = 0.16) +
  geom_smooth(data = d_f1[Subtype_Group == "TNBC"], method = "lm", se = FALSE, color = pal_cns[["red"]], linewidth = 0.90) +
  geom_point(size = 2.4, alpha = 0.88) +
  annotate("text", x = max(d_f1$CXCL8, na.rm=TRUE)*0.98, y = max(d_f1$LN_IC50, na.rm=TRUE)*0.98,
           label = f1_lbl, hjust = 1, vjust = 1, size = 2.55, fontface = "bold", color = "#1E293B") +
  scale_color_manual(values = c("ER+" = pal_cns[["blue"]], "HER2+" = pal_cns[["amber"]], "TNBC" = pal_cns[["red"]])) +
  labs(title = "F1. CXCL8 and dasatinib",
       x = expression(bold(CXCL8~(log[2]~TPM + 1))), y = expression(bold(Dasatinib~ln~IC[50])), color = NULL) +
  theme_cns(8.8) + theme(legend.position = "bottom")

d_f2 <- gb_mat[DRUG_ID == 1200L]
f2_lbl <- spearman_label_26b("PPARG", 1200L, d_f2, "ER+ lines")

p8f2 <- ggplot(d_f2, aes(x = PPARG, y = LN_IC50, color = Subtype)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75, alpha = 0.16) +
  geom_smooth(data = d_f2[Subtype_Group == "ER_pos"], method = "lm", se = FALSE, color = pal_cns[["blue"]], linewidth = 0.90) +
  geom_point(size = 2.4, alpha = 0.88) +
  annotate("text", x = min(d_f2$PPARG, na.rm=TRUE) + 0.15, y = max(d_f2$LN_IC50, na.rm=TRUE)*0.98,
           label = f2_lbl, hjust = 0, vjust = 1, size = 2.55, fontface = "bold", color = "#1E293B") +
  scale_color_manual(values = c("ER+" = pal_cns[["blue"]], "HER2+" = pal_cns[["amber"]], "TNBC" = pal_cns[["red"]])) +
  labs(title = "F2. PPARG and fulvestrant",
       x = expression(bold(PPARG~(log[2]~TPM + 1))), y = expression(bold(Fulvestrant~ln~IC[50])), color = NULL) +
  theme_cns(8.8) + theme(legend.position = "bottom")

fig8_r2 <- patchwork::wrap_elements(full = (p8b | km_8c) + plot_layout(widths = c(1.08, 1)))
fig8_r3 <- patchwork::wrap_elements(full = (p8d1 | p8d2))
fig8_r5 <- patchwork::wrap_elements(full = (p8f1 | p8f2))
fig8_all <- (p8a / fig8_r2 / fig8_r3 / p8e / fig8_r5) +
  plot_layout(heights = c(1.03, 1.30, 1.31, 1.42, 1.22))
save_frontiers_26b(fig8_all, "Fig8_Therapy_Response_ISPY2_and_GDSC2", height_mm = 253)
rm(fig8_all, fig8_r2, fig8_r3, fig8_r5, p8a, p8b, km_8c, p8d1, p8d2, p8e, p8f1, p8f2, gb_mat); gc()

# ==============================================================================
# FIG S1: 补全 A1-A6 / B1-B6 编号 + 4 方法共享图例 + 有限域拟合线 + 斑马纹 LOO
# ==============================================================================
cat("=== Building Fig S1 (MR Sensitivity with Full Legend & Subpanel Labels) ===\n")

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

mr_cols <- c("Inverse variance weighted" = pal_cns[["red"]], "Weighted median" = pal_cns[["green"]],
             "MR Egger" = pal_cns[["blue"]], "Weighted mode" = pal_cns[["amber"]])

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
    geom_hline(yintercept = 0, color = "#CBD5E1", linetype = "dotted") +
    geom_errorbar(aes(ymin = by - 1.96 * se.outcome, ymax = by + 1.96 * se.outcome), color = "#94A3B8", width = 0, linewidth = 0.4) +
    geom_errorbar(aes(xmin = bx - 1.96 * se.exposure, xmax = bx + 1.96 * se.exposure), orientation = "y", color = "#94A3B8", width = 0, linewidth = 0.4) +
    geom_segment(data = ms, aes(x = x1, xend = x2, y = y1, yend = y2, color = method), linewidth = 0.92) +
    geom_point(size = 2.1, color = "#0F172A") +
    scale_color_manual(values = mr_cols, breaks = names(mr_cols), name = "MR Method", drop = FALSE) +
    labs(title = sprintf("%s. %s", tk$tag, panel_label), x = "SNP effect on exposure", y = "SNP effect on BRCA") +
    theme_cns(8.9) + theme(legend.position = "bottom", plot.title.position = "plot",
                           plot.title = element_text(size = 9.3, lineheight = 1.00))
  sc_plots[[tk$lbl]] <- p_sc

  ld <- copy(loo_dt[id.exposure == tk$eid & id.outcome == tk$oid])
  ld[, is_all := (SNP == "All")]
  ld <- ld[order(is_all, -b)]
  ld[, snp_ord := factor(SNP, levels = unique(rev(SNP)))]
  ld[, y_num := as.numeric(snp_ord)]
  z_lo <- ld[y_num %% 2 == 1]

  p_lo <- ggplot(ld, aes(x = b, y = snp_ord, color = is_all)) +
    geom_rect(data = z_lo, aes(ymin = y_num - 0.46, ymax = y_num + 0.46, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
    geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), orientation = "y", width = 0.25, linewidth = 0.58) +
    geom_point(size = 2.1) +
    scale_color_manual(values = c("FALSE" = "#475569", "TRUE" = pal_cns[["red"]]), guide = "none") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
    labs(title = sprintf("%s. %s", tk$ltag, panel_label), x = expression("IVW " * beta), y = NULL) +
    theme_cns(8.9) + theme(plot.title.position = "plot", plot.title = element_text(size = 9.3, lineheight = 1.00))
  loo_plots[[tk$lbl]] <- p_lo
}

s1_leg <- cowplot::get_legend(sc_plots[[1]] + theme(legend.position = "bottom",
                                                    legend.title = element_text(face = "bold", size = 9.5),
                                                    legend.text = element_text(size = 9.0)))
sc_clean <- lapply(sc_plots, function(p) p + theme(legend.position = "none"))
s1_grid  <- wrap_plots(c(sc_clean, loo_plots), ncol = 3)
fig_s1   <- s1_grid / patchwork::wrap_elements(full = s1_leg) + plot_layout(heights = c(1.0, 0.055)) +
  plot_annotation(
    title = "Supplementary Figure S1. MR sensitivity and leave-one-out analyses",
    theme = theme(plot.title = element_text(face = "bold", size = 12.2, family = "sans", color = "#0F172A"))
  )
save_frontiers_26b(fig_s1, "FigS1_Upstream_MR_Sensitivity", height_mm = 193)
rm(fig_s1, sc_plots, sc_clean, loo_plots, s1_grid); gc()

# ==============================================================================
# FIG S2: 斑马纹异质性分解 (S2A) + 跨队列欧洲菌群 Meta (S2B)
# ==============================================================================
cat("=== Building Fig S2 (Zebra Heterogeneity + Microbiome Meta) ===\n")

het_map <- c(
  "All_5_Cohorts_Overall_OS"        = "All 5 Cohorts (OS)",
  "LOO_Excluding_METABRIC"          = "Excl. METABRIC",
  "LOO_Excluding_GSE42568"          = "Excl. GSE42568",
  "LOO_Excluding_TCGA_BRCA"         = "Excl. TCGA-BRCA",
  "LOO_Excluding_GSE20685"          = "Excl. GSE20685",
  "LOO_Excluding_GSE7390"           = "Excl. GSE7390",
  "Microarray_4Cohorts_Overall_OS"  = "Microarray Cohorts Only (n = 4)",
  "GEO_3Cohorts_Affy_OS"            = "Affymetrix GEO (k = 3)",
  "GEO_Large_GSE20685_GSE7390_OS"   = "Large Affymetrix GEO (GSE20685 + GSE7390)",
  "ER_neg_Stratum_TCGA_METABRIC_OS" = "TCGA + METABRIC ER\u2212 Stratum"
)

het_dt <- fread(file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))
het_dt[, `:=`(row_label = ifelse(analysis %in% names(het_map), unname(het_map[analysis]), gsub("_", " ", analysis)),
              effect = RE_HR, low = RE_LCI, high = RE_UCI,
              stats = sprintf("%.2f [%.2f, %.2f] | %.3g | %.0f%%", RE_HR, RE_LCI, RE_UCI, RE_pval, I2_pct))]
ps2a <- forest_26b(het_dt, "A. CXCL8 overall survival: heterogeneity and sensitivity", "Random-effects HR (log scale)",
                   1, "HR [95% CI] | P | I\u00b2", c(1.7, 1.1, 1.15), log_x = TRUE)

# 同秩主表中的所有配对。随机效应 OR、CI、P 和 BH 均取自相同行。
cc_meta <- fread(file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))[analysis_scope == "same_rank"]
if (!nrow(cc_meta) || anyDuplicated(cc_meta, by = c("core_taxon", "outcome_label")) ||
    !all(is.finite(unlist(cc_meta[, .(RE_OR, RE_LCI, RE_UCI, RE_pval, RE_FDR)])))) {
  stop("Missing or ambiguous same-rank random-effects microbiome results")
}
cc_meta[, outcome_short := clean_subtype_str(outcome_label)]
cc_meta[, outcome_order := match(outcome_short, c("Overall", "ER+", "ER\u2212"))]
setorder(cc_meta, core_taxon, outcome_order)
cc_meta[, `:=`(row_label = sprintf("%s | %s | k=%d", core_taxon, outcome_short, k_cohorts),
               effect = RE_OR, low = RE_LCI, high = RE_UCI,
               stats = sprintf("%.3f [%.3f, %.3f] | %.3g | %.3f | %.0f%%", RE_OR, RE_LCI, RE_UCI, RE_pval, RE_FDR, I2_pct))]
ps2b <- forest_26b(cc_meta, "B. Cross-cohort microbiome MR meta-analysis",
                   "Random-effects OR (log scale)", 1, "OR [95% CI] | P | BH q | I\u00b2",
                   c(1.43, 1.01, 1.51), log_x = TRUE,
                   subtitle = "Same taxonomic rank. Exploratory estimates (k = number of cohorts).")
fig_s2 <- ps2a / ps2b + plot_layout(heights = c(1, 1.65)) +
  plot_annotation(title = "Supplementary Figure S2. Heterogeneity and microbiome meta-analysis",
                  theme = theme(plot.title = element_text(face = "bold", size = 11.1, color = "#0F172A")))
save_frontiers_26b(fig_s2, "FigS2_HeterogeneityLOO_and_MicrobiomeMeta", height_mm = 182)
rm(fig_s2, ps2a, ps2b); gc()

# ==============================================================================
# FIG S3: 极简角标外部双单细胞验证 + HPA 参考单细胞(含数值) + 95% CI GSE9195 KM
# ==============================================================================
cat("=== Building Fig S3 (Minimalist External scRNA + HPA nCPM + 95% CI KM) ===\n")

# 只复用 26a 的两个函数定义，避免重画其他脚本负责的 Fig 6/7。
umap_helpers <- c("plot_umap_cns_ct", "plot_umap_cns_gene")
umap_env <- environment()
for (expr in parse(file.path(plot_script_dir, "26a_generate_cns_figs6to7.R"), encoding = "UTF-8")) {
  if (is.call(expr) && identical(expr[[1L]], quote(`<-`)) && is.symbol(expr[[2L]]) &&
      as.character(expr[[2L]]) %in% umap_helpers && is.call(expr[[3L]]) &&
      identical(expr[[3L]][[1L]], quote(`function`))) eval(expr, envir = umap_env)
}
stopifnot(all(vapply(umap_helpers, exists, logical(1), envir = umap_env, mode = "function", inherits = FALSE)))
sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
em_df <- copy(sc_cache[["BRCA_EMTAB8107"]])
ct_em_clean <- clean_lineage_str(em_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_em <- setNames(colorRampPalette(c(pal_cns[["blue"]], pal_cns[["red"]], pal_cns[["green"]], pal_cns[["cyan"]],
                                      pal_cns[["coral"]], pal_cns[["slate"]], pal_cns[["amber"]], "#7C3AED", "#EC4899", "#14B8A6", "#64748B"))(length(ct_em_clean)), ct_em_clean)
ps3a1 <- plot_umap_cns_ct(em_df, pal_em, "A1. E-MTAB-8107 Lineages (n = 33,043)")
ps3a2 <- plot_umap_cns_gene(em_df, "PPARG", "A2. E-MTAB-8107: PPARG", pal_mid = "#60A5FA", pal_hi = "#1E3A8A")
ps3a3 <- plot_umap_cns_gene(em_df, "CXCL8", "A3. E-MTAB-8107: CXCL8", pal_mid = "#F59E0B", pal_hi = "#9F1239")

g14_df <- copy(sc_cache[["BRCA_GSE148673"]])
ct_g14_clean <- clean_lineage_str(g14_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_g14 <- setNames(colorRampPalette(c(pal_cns[["blue"]], pal_cns[["red"]], pal_cns[["green"]], pal_cns[["cyan"]],
                                       pal_cns[["coral"]], pal_cns[["slate"]], pal_cns[["amber"]], "#7C3AED", "#14B8A6"))(length(ct_g14_clean)), ct_g14_clean)
ps3b1 <- plot_umap_cns_ct(g14_df, pal_g14, "B1. GSE148673 Lineages (n = 10,359)")
ps3b2 <- plot_umap_cns_gene(g14_df, "PPARG", "B2. GSE148673: PPARG", pal_mid = "#60A5FA", pal_hi = "#1E3A8A")
ps3b3 <- plot_umap_cns_gene(g14_df, "IL1B",  "B3. GSE148673: IL1B",  pal_mid = "#F59E0B", pal_hi = "#9F1239")

sc_hpa <- fread(file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
sc_keep <- c("Neutrophils","Macrophages","Fibroblasts","Adipocytes","Lymphatic endothelial cells","Vascular endothelial cells","Breast secretory cells","T-cells","B-cells")
sc_sub <- sc_hpa[gene_symbol %in% c("CXCL8","PPARG") & cell_type %in% sc_keep]
sc_sub[, cell_clean := tools::toTitleCase(gsub("-", " ", cell_type))]
sc_sub[, log10_ncpm := log10(nCPM + 1)]
sc_ord <- sc_sub[, .(mx = max(log10_ncpm)), by = cell_clean][order(mx)]$cell_clean
sc_sub[, cell_clean := factor(cell_clean, levels = sc_ord)]
sc_sub[, y_num := as.numeric(cell_clean)]
z_s3c <- data.table(y_num = seq_along(sc_ord))[y_num %% 2 == 1]

ps3c <- ggplot(sc_sub, aes(x = log10_ncpm, y = cell_clean, fill = gene_symbol)) +
  geom_rect(data = z_s3c, aes(ymin = y_num - 0.46, ymax = y_num + 0.46, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66) +
  geom_text(aes(label = ifelse(nCPM > 0, format(round(nCPM, 0), big.mark = ","), "0"),
                x = log10_ncpm + 0.05),
            position = position_dodge(width = 0.74), hjust = 0, size = 2.15, color = "#334155") +
  scale_fill_manual(values = c("CXCL8" = pal_cns[["red"]], "PPARG" = pal_cns[["blue"]])) +
  coord_cartesian(xlim = c(0, 4.35)) +
  labs(title = "C. HPA single-cell expression", x = expression(bold(log[10](nCPM + 1))), y = NULL, fill = NULL) +
  theme_cns(9.0) + theme(legend.position = "bottom")

g91_km <- copy(th_pc$gse9195)[!is.na(t_dmfs) & t_dmfs > 0 & !is.na(e_dmfs)]
g91_km[, Group := ifelse(PPARG >= median(PPARG), "High", "Low")]
km_s3d <- plot_km_cns(
  g91_km, "t_dmfs", "e_dmfs", "Group",
  "D. GSE9195 DMFS (ER+)", "Years", "DMFS probability", "PPARG",
  c("Low", "High"), c(pal_cns[["red"]], pal_cns[["blue"]]),
  pval_coord = c(0.4, 0.16), break_by = 3, max_x = 12, base_size = 9.0
)

figs3_r1 <- (ps3a1 | ps3a2 | ps3a3)
figs3_r2 <- (ps3b1 | ps3b2 | ps3b3)
figs3_r3 <- patchwork::wrap_elements(full = (ps3c | km_s3d) + plot_layout(widths = c(1.05, 1.05)))
fig_s3 <- figs3_r1 / figs3_r2 / figs3_r3 + plot_layout(heights = c(0.95, 0.95, 1.15)) +
  plot_annotation(
    title = "Supplementary Figure S3. Single-cell expression and tamoxifen outcomes",
    theme = theme(plot.title = element_text(face = "bold", size = 12.0, family = "sans", color = "#0F172A"))
  )
save_frontiers_26b(fig_s3, "FigS3_scRNA_ExternalCohorts_and_GSE9195", height_mm = 184)

cat("\n26b: Fig 8 and Fig S1–S3 exported.\n")
