# 02_fig5_bulk_survival.R
# Figure 5 (was Old Fig 4): Bulk tumor expression (Tumor vs Normal, PAM50),
# representative KM survival curves (fixed P/HR & risk-table spacing), and standardized OS meta-analysis.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Final Fig 5 (Bulk Expression & Survival Meta) ===\n")

tc      <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb      <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20     <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42     <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc   <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
diff_dt <- fread(file.path(res_dir, "module4_differential_expression.csv"))
s13     <- fread(file.path(res_dir, "planA_pam50_subtypes_stats.csv"))

# 5A: Tumor vs Normal PPARG
df_tn <- rbind(
  data.table(Cohort = "GSE42568",
             G_simple = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal", "Tumor"),
             Expr = as.numeric(g42$expr_target["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Normal", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Tumor",  Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)
df_tn[, G_simple := factor(G_simple, levels = c("Normal", "Tumor"))]
vio_5a <- build_half_violin_dt(df_tn, "G_simple", "Expr", facet_col = "Cohort", width = 0.30, side = "right")

p5a_gse <- diff_dt[cohort == "GSE42568"  & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
p5a_tcg <- diff_dt[cohort == "TCGA_BRCA" & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
brk_5a <- df_tn[, .(ymax = max(Expr, na.rm = TRUE), ymin = min(Expr, na.rm = TRUE)), by = Cohort]
brk_5a[, `:=`(y_brk = ymax + 0.08 * (ymax - ymin), y_tick = ymax + 0.03 * (ymax - ymin),
              y_txt = ymax + 0.14 * (ymax - ymin),
              p_lbl = c(paste0("P = ", format_p_clean(p5a_gse)), paste0("P = ", format_p_clean(p5a_tcg))))]

p5a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple, color = G_simple)) +
  geom_polygon(data = vio_5a, aes(x = vx, y = vy, group = grp_id, fill = G_simple),
               inherit.aes = FALSE, alpha = 0.32, color = pal_rev$slate, linewidth = 0.30) +
  ggrastr::rasterise(geom_jitter(aes(x = as.numeric(G_simple) - 0.16), width = 0.06, size = 0.45, alpha = 0.28), dpi = 300) +
  geom_boxplot(aes(x = as.numeric(G_simple) - 0.16), width = 0.14, outlier.shape = NA, alpha = 0.85, color = pal_rev$ink, linewidth = 0.40) +
  geom_segment(data = brk_5a, aes(x = 1, xend = 2, y = y_brk, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_segment(data = brk_5a, aes(x = 1, xend = 1, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_segment(data = brk_5a, aes(x = 2, xend = 2, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_text(data = brk_5a, aes(x = 1.5, y = y_txt, label = p_lbl), inherit.aes = FALSE, size = 2.35, color = pal_rev$ink) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_x_discrete(limits = c("Normal", "Tumor")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  scale_fill_manual(values = c("Normal" = "#38BDF8", "Tumor" = pal_rev$cxcl8)) +
  scale_color_manual(values = c("Normal" = "#0284C7", "Tumor" = pal_rev$cxcl8)) +
  labs(title = "A   PPARG (Tumor vs. Normal)", x = NULL, y = "Normalized log2 expr.") +
  theme_rev(8.2) + theme(legend.position = "none")

# 5B: PAM50 Subtypes (TCGA-BRCA)
pam_tc <- rbind(
  data.table(Gene = "CXCL8", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ])),
  data.table(Gene = "PPARG", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)[Subtype %in% c("Normal", "LumA", "LumB", "Her2", "Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"),
                           labels = c("Norm.", "LumA", "LumB", "HER2+", "Basal"))]
pam_tc[, Gene := factor(Gene, levels = c("CXCL8", "PPARG"))]
vio_5b <- build_half_violin_dt(pam_tc, "Subtype", "Expr", facet_col = "Gene", width = 0.32, side = "right")

kw_cx <- s13[cohort == "TCGA_BRCA" & gene_symbol == "CXCL8"]$kw_pval[1]
w_cx  <- wilcox.test(pam_tc[Gene == "CXCL8" & Subtype == "Basal"]$Expr, pam_tc[Gene == "CXCL8" & Subtype == "LumA"]$Expr)$p.value
kw_pp <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$kw_pval[1]
w_pp  <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$wilcox_Basal_vs_LumA[1]

ann_5b <- data.table(
  Gene = factor(c("CXCL8", "PPARG"), levels = c("CXCL8", "PPARG")),
  ymax = c(max(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), max(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  ymin = c(min(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), min(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  kw_lbl = c(paste0("KW P = ", format_p_clean(kw_cx)), paste0("KW P = ", format_p_clean(kw_pp))),
  w_lbl  = c(paste0("Basal vs LumA: P ", ifelse(w_cx < 2.2e-16, "< 2.2e-16", paste0("= ", format_p_clean(w_cx)))),
             paste0("Basal vs LumA: P = ", format_p_clean(w_pp)))
)
ann_5b[, `:=`(y_brk = ymax + 0.08 * (ymax - ymin), y_tick = ymax + 0.03 * (ymax - ymin),
              y_w = ymax + 0.16 * (ymax - ymin), y_kw = ymax + 0.30 * (ymax - ymin))]

p5b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_polygon(data = vio_5b, aes(x = vx, y = vy, group = grp_id, fill = Subtype),
               inherit.aes = FALSE, alpha = 0.32, color = pal_rev$slate, linewidth = 0.28) +
  geom_boxplot(aes(x = as.numeric(Subtype) - 0.14), width = 0.20, outlier.size = 0.35, alpha = 0.85, color = pal_rev$ink, linewidth = 0.38) +
  geom_segment(data = ann_5b, aes(x = 2, xend = 5, y = y_brk, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_segment(data = ann_5b, aes(x = 2, xend = 2, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_segment(data = ann_5b, aes(x = 5, xend = 5, y = y_tick, yend = y_brk), inherit.aes = FALSE, color = pal_rev$ink, linewidth = 0.42) +
  geom_text(data = ann_5b, aes(x = 3.5, y = y_w, label = w_lbl), inherit.aes = FALSE, size = 2.05, color = pal_rev$ink) +
  geom_text(data = ann_5b, aes(x = 3.0, y = y_kw, label = kw_lbl), inherit.aes = FALSE, size = 2.15, fontface = "bold", color = pal_rev$subink) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_x_discrete(limits = c("Norm.", "LumA", "LumB", "HER2+", "Basal")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  scale_fill_manual(values = c("Norm." = "#94A3B8", "LumA" = pal_rev$pparg, "LumB" = "#6366F1", "HER2+" = pal_rev$amber, "Basal" = pal_rev$cxcl8)) +
  labs(title = "B   PAM50 subtypes (TCGA-BRCA)", x = NULL, y = "log2(RSEM + 1)") +
  theme_rev(8.2) +
  theme(legend.position = "none", strip.text = element_text(face = "bold.italic"),
        axis.text.x = element_text(size = 6.9, angle = 22, hjust = 1))

# 5C–5F: Clean KM curves (zero P/HR overlap, zero risk-table overlap)
idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp],
                      expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5c <- plot_km_rev(df_tc_p, "time", "status", "Group", "C   TCGA ER+ OS (PPARG)", "Months", "Overall survival",
                     "PPARG", c("Low", "High"), c(pal_rev$cxcl8, pal_rev$pparg), break_by = 50, max_x = 250)

idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp],
                      expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5d <- plot_km_rev(df_mb_p, "time", "status", "Group", "D   METABRIC ER+ OS (PPARG)", "Months", "Overall survival",
                     "PPARG", c("Low", "High"), c(pal_rev$cxcl8, pal_rev$pparg), break_by = 100, max_x = 300)

pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5e <- plot_km_rev(df_g20_c, "time", "status", "Group", "E   GSE20685 OS (CXCL8)", "Years", "Overall survival",
                     "CXCL8", c("Low", "High"), c(pal_rev$pparg, pal_rev$cxcl8), break_by = 3, max_x = 15)

idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern],
                      expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5f <- plot_km_rev(df_mb_c, "time", "status", "Group", "F   METABRIC ER\u2212 OS (CXCL8)", "Months", "Overall survival",
                     "CXCL8", c("Low", "High"), c(pal_rev$pparg, pal_rev$cxcl8), break_by = 100, max_x = 300)
rm(tc, mb, g20, g42, cl_tc); gc()

# 5G: Standardized OS random-effects meta-analysis (Table-aligned forest)
meta_hr <- fread(file.path(res_dir, "planA_standardized_os_meta.csv"))
hr_sub <- meta_hr[(gene_symbol == "CXCL8" & stratum %chin% c("Overall", "ER_neg")) |
                  (gene_symbol == "PPARG" & stratum %chin% c("Overall", "ER_pos"))]
hr_sub[, `:=`(
  row_label = sprintf("%s | %s (%d cohorts, n=%s)", gene_symbol, clean_subtype_rev(stratum), k, format(total_n, big.mark = ",", trim = TRUE)),
  effect = RE_HR, low = RE_LCI, high = RE_UCI,
  col_group = ifelse(gene_symbol == "CXCL8", pal_rev$cxcl8, pal_rev$pparg),
  stats = sprintf("%.3f [%.3f, %.3f] | P=%.3g | I\u00b2=%.0f%%", RE_HR, RE_LCI, RE_UCI, RE_pval, I2_pct)
)]
p5g <- table_forest_rev(hr_sub, "G   Standardized overall survival meta-analysis (random effects)",
                        "Hazard ratio per SD (95% Wald CI)", 1, "HR [95% CI] | P | I\u00b2",
                        widths = c(1.65, 1.30, 1.65), log_x = TRUE,
                        x_range = c(0.72, 1.50), x_breaks = c(0.8, 1.0, 1.2, 1.4))

fig5_r1 <- wrap_elements(full = (p5a | p5b) + plot_layout(widths = c(0.95, 1.15)))
fig5_all <- wrap_plots(fig5_r1, wrap_elements(full = (km_5c | km_5d)),
                       wrap_elements(full = (km_5e | km_5f)), p5g,
                       ncol = 1, heights = c(56, 63, 63, 46))
save_rev_figure(fig5_all, "Fig5_Bulk_Expression_and_Survival_Meta", width_mm = 180, height_mm = 230, dpi_png = 350)
cat("Fig 5 complete.\n")
