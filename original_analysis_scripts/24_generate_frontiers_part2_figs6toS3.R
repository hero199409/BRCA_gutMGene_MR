# 24_generate_frontiers_part2_figs6toS3.R
# Nature / Cell / Frontiers 顶刊美学重塑版 Part 2 (Fig 6 - Fig 8 + Fig S1 - Fig S3)
# 特点：零副标题小字、零代码下划线/ASCII符号、开放式L型主题、无框圆形HPA切片配底注、I-SPY2显著性括号、严格语义着色

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/23_0_plot_helpers.R"))

sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
sc_sum   <- fread(file.path(res_dir, "planA_scrna_celltype_summary.csv"))

plot_umap_ct <- function(umap_dt, ct_pal, ttl, pt_sz = 0.22) {
  d <- copy(umap_dt)
  d[, Lineage_clean := clean_lineage_str(Celltype_major)]
  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = Lineage_clean)) +
    ggrastr::rasterise(geom_point(size = pt_sz, alpha = 0.65), dpi = 320) +
    scale_color_manual(values = ct_pal, guide = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    labs(title = ttl, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    theme_journal(8.8) + theme(legend.position = "right", legend.key.size = unit(0.32, "cm"))
}

plot_umap_gene <- function(umap_dt, gene, ttl, hi_col = "#DC2626", pt_sz = 0.22) {
  d <- copy(umap_dt[, .(UMAP_1, UMAP_2, val = get(gene))])[order(val)]
  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = val)) +
    ggrastr::rasterise(geom_point(size = pt_sz, alpha = ifelse(d$val > 0, 0.85, 0.25)), dpi = 320) +
    scale_color_gradient(low = "#E2E8F0", high = hi_col) +
    labs(title = ttl, x = "UMAP 1", y = "UMAP 2", color = "Expr") +
    theme_journal(8.8) + theme(legend.position = "right", legend.key.height = unit(0.55, "cm"))
}

# ==============================================================================
# FIG 6: 肿瘤单细胞定位 (GSE176078) + 亚型解卷积 + 去边框大尺寸 HPA 免疫组化切片
# ==============================================================================
sc_df    <- copy(sc_cache[["BRCA_GSE176078"]])
sc_umap  <- sc_df[!is.na(UMAP_1) & !is.na(UMAP_2)]
sc_sum_p <- copy(sc_sum[dataset == "BRCA_GSE176078"])

ct_raw6 <- sc_umap[, .N, by = Celltype_major][order(-N)]$Celltype_major
ct_clean6 <- clean_lineage_str(ct_raw6)
ct_pal6 <- setNames(colorRampPalette(c(pal_npg[["blue"]], pal_npg[["red"]], pal_npg[["green"]], pal_npg[["cyan"]],
                                       pal_npg[["coral"]], pal_npg[["slate"]], pal_npg[["amber"]], "#7C3AED", "#EC4899"))(length(ct_clean6)), ct_clean6)

p6a_ct <- plot_umap_ct(sc_umap, ct_pal6, "A1. Tumor scRNA Lineages (GSE176078)")
p6a_cx <- plot_umap_gene(sc_umap, "CXCL8", "A2. CXCL8 Expression", hi_col = pal_npg[["red"]])
p6a_pp <- plot_umap_gene(sc_umap, "PPARG", "A3. PPARG Expression", hi_col = pal_npg[["blue"]])
p6a_pl <- plot_umap_gene(sc_umap, "PLIN1", "A4. Adipocyte Control (PLIN1)", hi_col = "#64748B")

dot_genes6 <- c("CXCL8","IL1B","FFAR2","PPARG","FABP4","ADIPOQ","PLIN1","LEP","EPCAM","VWF","CD68")
dot_dt6 <- copy(sc_sum_p[gene %in% dot_genes6])
dot_dt6[, gene := factor(gene, levels = rev(dot_genes6))]
ct_ord6 <- clean_lineage_str(c("Malignant","Epithelial","Endothelial","PVL","Fibroblasts","Mono/Macro","T_cells","B","Plasma"))
dot_dt6[, Lineage := factor(clean_lineage_str(celltype_major), levels = ct_ord6)]

p6b <- ggplot(dot_dt6, aes(x = Lineage, y = gene, size = pct_expr, color = mean_expr)) +
  geom_point() +
  scale_size(range = c(0.5, 6.2), name = "% Expressing") +
  scale_color_gradient(low = "#E0F2FE", high = "#BE123C", name = "Mean Expr") +
  labs(title = "B. Lineage Specificity DotPlot (GSE176078)", x = NULL, y = NULL) +
  theme_panel_grid(9.2) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

st_cmp <- fread(file.path(res_dir, "planA_scrna_subtype_comparison.csv"))
st_sub <- st_cmp[celltype %in% c("ALL_CELLS", "Endothelial", "Malignant", "Mono/Macro") & gene %in% c("PPARG", "CXCL8")]
st_sub[, celltype_lbl := factor(celltype, levels = c("ALL_CELLS", "Endothelial", "Malignant", "Mono/Macro"),
                                labels = c("All Cells", "Endothelial", "Malignant Cells", "Mono/Macrophages"))]
st_sub[, subtype := factor(subtype, levels = c("ER+", "HER2+", "TNBC"))]
st_sub[, gene := factor(gene, levels = c("PPARG", "CXCL8"))]

p6c <- ggplot(st_sub, aes(x = subtype, y = mean_expr, fill = subtype)) +
  geom_col(width = 0.60) +
  geom_text(aes(label = sprintf("%.2f", mean_expr)), vjust = -0.25, size = 2.35) +
  facet_grid(gene ~ celltype_lbl, scales = "free_y") +
  scale_fill_manual(values = c("ER+" = pal_npg[["blue"]], "HER2+" = pal_npg[["amber"]], "TNBC" = pal_npg[["red"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "C. Molecular Subtype Cellular Deconvolution", x = NULL, y = "Mean Normalized Expression") +
  theme_journal(9.0) + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

# 无灰边白框的圆形 HPA IHC 病理切片 + 底部简洁图注
p6d1 <- img_panel(file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"),
                  bottom_caption = "D1. Normal Breast \u2014 PPARG (Adipocytes Low; Ductal ND)", show_border = FALSE)
p6d2 <- img_panel(file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"),
                  bottom_caption = "D2. Breast Carcinoma \u2014 PPARG (4 Med / 6 Low / 2 ND)", show_border = FALSE)
p6d3 <- img_panel(file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"),
                  bottom_caption = "D3. Breast Carcinoma \u2014 CXCL8 (12/12 ND, Secreted)", show_border = FALSE)

fig6_r1 <- (p6a_ct | p6a_cx | p6a_pp | p6a_pl)
fig6_r2 <- (p6b | p6c) + plot_layout(widths = c(1.0, 1.15))
fig6_r3 <- (p6d1 | p6d2 | p6d3)
fig6_all <- fig6_r1 / fig6_r2 / fig6_r3 + plot_layout(heights = c(0.95, 1.05, 1.05))
save_pub(fig6_all, "Fig6_Tumor_scRNA_and_HPA_Pathology", w = 17.5, h = 15.8)
rm(fig6_all, fig6_r1, fig6_r2, fig6_r3, p6a_ct, p6a_cx, p6a_pp, p6a_pl, p6b, p6c, p6d1, p6d2, p6d3, sc_df, sc_umap); gc()

# ==============================================================================
# FIG 7: TME 免疫微环境 (洁净英文特征名 + ER+) + Title Case GSEA + 3队列跨库一致性
# ==============================================================================
tme_map <- c(
  "quanTIseq__B cell"                = "B Cells (quanTIseq)",
  "quanTIseq__Macrophage M1"         = "M1 Macrophages (quanTIseq)",
  "MCPcounter__Neutrophil"           = "Neutrophils (MCPcounter)",
  "quanTIseq__NK cell"               = "NK Cells (quanTIseq)",
  "TMB_NONSYNONYMOUS"                = "Nonsynonymous TMB",
  "MCPcounter__T cell CD8+"          = "CD8+ T Cells (MCPcounter)",
  "quanTIseq__Neutrophil"            = "Neutrophils (quanTIseq)",
  "quanTIseq__T cell CD8+"           = "CD8+ T Cells (quanTIseq)",
  "MCPcounter__T cell"               = "T Cells (MCPcounter)",
  "MCPcounter__B cell"               = "B Cells (MCPcounter)",
  "MCPcounter__Endothelial cell"     = "Endothelial Cells (MCPcounter)",
  "ESTIMATE__stroma score"           = "Stromal Score (ESTIMATE)",
  "ESTIMATE__tumor purity"           = "Tumor Purity (ESTIMATE)",
  "ESTIMATE__estimate score"         = "ESTIMATE Score",
  "MCPcounter__cytotoxicity score"   = "Cytotoxicity Score (MCPcounter)",
  "quanTIseq__uncharacterized cell"  = "Uncharacterized Cells (quanTIseq)",
  "MCPcounter__Myeloid dendritic cell" = "Myeloid Dendritic Cells (MCPcounter)"
)

cor5 <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
tme_top <- cor5[((gene_symbol=="CXCL8" & stratum=="ER_neg") | (gene_symbol=="PPARG" & stratum=="ER_pos")) &
                  feature_type %in% c("TME_Cell_Score","Checkpoint","Genomic_Feature") & FDR < 0.05][
                    order(gene_symbol, -abs(rho))][, .SD[1:min(10, .N)], by = gene_symbol]
tme_top[, feat_clean := ifelse(feature %in% names(tme_map), unname(tme_map[feature]), gsub("_", " ", feature))]
tme_top[, facet_lbl := ifelse(gene_symbol=="CXCL8", "A1. CXCL8 in ER\u2212 BRCA (n = 237)", "A2. PPARG in ER+ BRCA (n = 783)")]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feat_clean)]
tme_top <- tme_top[order(facet_lbl, rho)]; tme_top[, feat_ord := factor(feat_ord, levels = unique(feat_ord))]

p7a <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_segment(aes(x = 0, xend = rho, y = feat_ord, yend = feat_ord), linewidth = 0.8) +
  geom_point(size = 2.75) + facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("Positive" = pal_npg[["red"]], "Negative" = pal_npg[["blue"]])) +
  coord_cartesian(xlim = c(-0.52, 0.66)) +
  labs(title = "A. Tumor Microenvironment & Genomic Correlations", x = "Spearman \u03c1", y = NULL, color = NULL) +
  theme_journal(9.3) + theme(legend.position = "bottom")

gsea_map <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"           = "TNF\u03b1\u2013NF-\u03baB Signaling",
  "HALLMARK_MTORC1_SIGNALING"                  = "mTORC1 Signaling",
  "HALLMARK_G2M_CHECKPOINT"                    = "G2/M Checkpoint",
  "HALLMARK_HYPOXIA"                           = "Hypoxia",
  "HALLMARK_GLYCOLYSIS"                        = "Glycolysis",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Inflammatory Response",
  "HALLMARK_E2F_TARGETS"                       = "E2F Targets",
  "HALLMARK_ANGIOGENESIS"                      = "Angiogenesis",
  "HALLMARK_ALLOGRAFT_REJECTION"               = "Allograft Rejection",
  "HALLMARK_KRAS_SIGNALING_UP"                 = "KRAS Signaling Up",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "IL-6/JAK/STAT3 Signaling",
  "HALLMARK_COAGULATION"                       = "Coagulation",
  "HALLMARK_MYOGENESIS"                        = "Myogenesis"
)

gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene=="CXCL8" & stratum=="ER_neg") | (target_gene=="PPARG" & stratum=="ER_pos")) & p.adjust < 0.05][
  order(target_gene, -abs(NES))][, .SD[1:min(9, .N)], by = target_gene]
gs_top[, clean_id := ifelse(ID %in% names(gsea_map), unname(gsea_map[ID]), gsub("_", " ", sub("^HALLMARK_", "", ID)))]
gs_top[, facet_lbl := ifelse(target_gene=="CXCL8", "B1. CXCL8 (ER\u2212 BRCA)", "B2. PPARG (ER+ BRCA)")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]; gs_top[, id_ord := factor(id_ord, levels = id_ord)]

p7b <- ggplot(gs_top, aes(x = NES, y = id_ord, fill = -log10(p.adjust))) +
  geom_col(width = 0.72) + facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_fill_gradient(low = "#93C5FD", high = "#1E3A8A") +
  labs(title = "B. Hallmark Pathway Enrichment (GSEA)",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = expression(bold(-log[10](FDR)))) +
  theme_journal(9.3) + theme(legend.position = "bottom")

shared_ct <- c("Endothelial", "Mono/Macro", "Fibroblasts", "Malignant")
rep_genes <- c("PPARG", "FABP4", "CXCL8", "IL1B", "PLIN1", "ADIPOQ")
grid_all  <- CJ(dataset = c("BRCA_GSE176078", "BRCA_EMTAB8107", "BRCA_GSE148673"),
                celltype_major = shared_ct, gene = rep_genes)
rep_dt <- merge(grid_all, sc_sum[, .(dataset, celltype_major, gene, pct_expr, mean_expr)],
                by = c("dataset", "celltype_major", "gene"), all.x = TRUE)
rep_dt[, ds_lbl := factor(dataset,
                          levels = c("BRCA_GSE176078", "BRCA_EMTAB8107", "BRCA_GSE148673"),
                          labels = c("GSE176078 (100,064 cells)", "E-MTAB-8107 (33,043 cells)", "GSE148673 (10,359 cells)"))]
rep_dt[, Lineage := factor(clean_lineage_str(celltype_major), levels = clean_lineage_str(shared_ct))]
rep_dt[, gene := factor(gene, levels = rep_genes)]
rep_dt[, is_na := is.na(pct_expr)]
rep_dt[, pct_lbl := ifelse(is_na, "ND", ifelse(pct_expr >= 1, sprintf("%.1f%%", pct_expr), sprintf("%.2f%%", pct_expr)))]

p7c <- ggplot(rep_dt, aes(x = Lineage, y = ds_lbl)) +
  geom_tile(data = rep_dt[is_na == TRUE], fill = "#F8FAFC", color = "#E2E8F0") +
  geom_point(data = rep_dt[is_na == FALSE], aes(size = pct_expr, color = mean_expr)) +
  geom_text(data = rep_dt[is_na == FALSE], aes(label = pct_lbl), position = position_nudge(y = 0.28), size = 2.4, fontface = "bold", color = "#334155") +
  geom_text(data = rep_dt[is_na == TRUE],  aes(label = pct_lbl), size = 2.4, fontface = "italic", color = "#94A3B8") +
  facet_wrap(~ gene, nrow = 1) +
  scale_size(range = c(1.0, 7.8), name = "% Expressing") +
  scale_color_gradient(low = "#BAE6FD", high = "#BE123C", name = "Mean Expr") +
  scale_y_discrete(limits = rev(c("GSE176078 (100,064 cells)", "E-MTAB-8107 (33,043 cells)", "GSE148673 (10,359 cells)")),
                   expand = expansion(mult = c(0.18, 0.34))) +
  labs(title = "C. Cross-Cohort Tumor scRNA Consistency (n = 143,466 Cells)", x = NULL, y = NULL) +
  theme_panel_grid(9.2) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "right",
        panel.spacing.x = unit(0.6, "lines"))

fig7_all <- (p7a | p7b) / p7c + plot_layout(heights = c(1.12, 0.95))
save_pub(fig7_all, "Fig7_TME_GSEA_and_CrossCohort_scRNA", w = 17.6, h = 13.5)
rm(fig7_all, p7a, p7b, p7c, rep_dt, grid_all); gc()

# ==============================================================================
# FIG 8: 真实临床治疗响应 (修复8A颜色一致性) + I-SPY2显式括号 + 修复GDSC2药敏
# ==============================================================================
th_dt <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))
th_pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

# 8A: 严格按所展示的 Cox P 值着色（修复 GSE9195 DMFS 颜色冲突），并显式注明 KM P=0.042
endo_dt <- copy(th_dt[gene == "PPARG" & cohort %in% c("METABRIC", "GSE9195", "GSE17705") & stratum == "ER_pos"])
endo_dt[, arm_clean := fifelse(treatment_arm == "Hormone_Therapy_YES", "Endocrine (+)",
                       fifelse(treatment_arm == "Hormone_Therapy_NO",  "Endocrine (\u2212)",
                       fifelse(treatment_arm == "Tamoxifen_5yr",       "Tamoxifen", treatment_arm)))]
endo_dt[, lbl := sprintf("%s | %s %s (n = %s)", cohort, arm_clean, endpoint, format(n_samples, big.mark = ","))]
endo_dt[, tag := ifelse(pval < 0.05, "Cox P < 0.05", "Cox P \u2265 0.05")]
endo_dt[, hr_anno := ifelse(cohort == "GSE9195" & endpoint == "DMFS",
                            sprintf("HR = %.3f (Cox P = %.3f; KM P = %.3f)", estimate, pval, km_or_wilcox_pval),
                            sprintf("HR = %.3f, P = %.3g", estimate, pval))]
endo_dt[, lbl := factor(lbl, levels = rev(lbl))]

p8a <- ggplot(endo_dt, aes(x = estimate, y = lbl, color = tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = lci95, xmax = uci95), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(size = 2.8) +
  geom_text(aes(x = 1.48, label = hr_anno), hjust = 0, size = 2.45, show.legend = FALSE) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1.0, 1.5)) +
  scale_color_manual(values = c("Cox P < 0.05" = pal_npg[["blue"]], "Cox P \u2265 0.05" = pal_npg[["null"]])) +
  coord_cartesian(xlim = c(0.08, 16.5)) +
  labs(title = "A. ER+ Endocrine Therapy Outcomes (PPARG)", x = "Hazard Ratio (95% CI, log scale)", y = NULL, color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

nact_pcr_dt <- copy(th_dt[gene == "CXCL8" & cohort %in% c("GSE25066", "GSE32646") & endpoint == "pCR_vs_RD" & stratum %in% c("Overall", "ER_neg", "TNBC")])
nact_pcr_dt[, lbl := sprintf("%s %s (n = %d, %d pCR)", cohort, clean_subtype_str(stratum), n_samples, n_events)]
nact_pcr_dt[, tag := ifelse(pval < 0.05, "P < 0.05 (Unadjusted)", "ER-Stratified (P \u2265 0.05)")]
nact_pcr_dt[, lbl := factor(lbl, levels = rev(lbl))]

p8b <- ggplot(nact_pcr_dt, aes(x = estimate, y = lbl, color = tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = lci95, xmax = uci95), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(size = 2.65) +
  geom_text(aes(x = uci95 + 0.03, label = sprintf("OR = %.2f, AUC = %.2f, P = %.3g", estimate, auc, pval)),
            hjust = 0, size = 2.45, show.legend = FALSE) +
  scale_color_manual(values = c("P < 0.05 (Unadjusted)" = pal_npg[["red"]], "ER-Stratified (P \u2265 0.05)" = "#64748B")) +
  coord_cartesian(xlim = c(0.55, 3.35)) +
  labs(title = "B. NACT Pathological Complete Response (CXCL8)", x = "Odds Ratio for pCR (95% CI)", y = NULL, color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

g25_km <- copy(th_pc$gse25066)[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs)]
g25_km[, Group := ifelse(CXCL8 >= median(CXCL8), "High", "Low")]
km_8c <- plot_km_with_risk(
  g25_km, "t_drfs", "e_drfs", "Group",
  "C. GSE25066 NACT DRFS (CXCL8)", "Years", "Distant Relapse-Free Survival", "CXCL8",
  c("Low", "High"), c(pal_npg[["blue"]], pal_npg[["red"]]),
  pval_coord = c(0.3, 0.16), break_by = 2, max_x = 8, base_size = 8.6
)

# 8D1 & 8D2: I-SPY2 免疫治疗 —— 箱线图上方加对比支架线并标注显式斜体 P 值
g173_df <- copy(th_pc$gse173839)[!is.na(pCR)]
g173_df[, Resp := factor(ifelse(pCR == 1, "pCR", "RD"), levels = c("RD", "pCR"))]
g173_df[, Arm_lbl := factor(ifelse(Arm == "durvalumab/olaparib", "Durvalumab + Olaparib (n = 71)", "Control Paclitaxel (n = 29)"),
                            levels = c("Control Paclitaxel (n = 29)", "Durvalumab + Olaparib (n = 71)"))]
y_d1_top <- max(g173_df$CXCL8, na.rm = TRUE)
d1_span  <- diff(range(g173_df$CXCL8, na.rm = TRUE))
brk_d1 <- data.table(
  Arm_lbl = factor(c("Control Paclitaxel (n = 29)", "Durvalumab + Olaparib (n = 71)"), levels = levels(g173_df$Arm_lbl)),
  y_top   = y_d1_top + 0.08 * d1_span,
  y_tick  = y_d1_top + 0.04 * d1_span,
  y_txt   = y_d1_top + 0.13 * d1_span,
  p_lbl   = c("italic(P) == 0.15", "italic(P) == 0.0029")
)

p8d1 <- ggplot(g173_df, aes(x = Resp, y = CXCL8, fill = Resp)) +
  geom_boxplot(width = 0.45, alpha = 0.85, outlier.size = 0.8) +
  geom_jitter(width = 0.13, size = 1.15, alpha = 0.55) +
  geom_segment(data = brk_d1, aes(x = 1, xend = 2, y = y_top, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_segment(data = brk_d1, aes(x = 1, xend = 1, y = y_tick, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_segment(data = brk_d1, aes(x = 2, xend = 2, y = y_tick, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_text(data = brk_d1, aes(x = 1.5, y = y_txt, label = p_lbl), inherit.aes = FALSE, parse = TRUE, size = 3.0, color = "#0F172A") +
  facet_wrap(~ Arm_lbl) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  scale_fill_manual(values = c("RD" = pal_npg[["null"]], "pCR" = pal_npg[["red"]])) +
  labs(title = "D1. I-SPY2 CXCL8 \u00d7 Immunotherapy Interaction (P_int = 0.0056)", x = NULL, y = "CXCL8 Expression") +
  theme_journal(9.0) + theme(legend.position = "none")
# 用正规 Unicode 下标样式修饰标题中的 P_int -> Interaction P = 0.0056
p8d1 <- p8d1 + labs(title = "D1. I-SPY2 CXCL8 \u00d7 Immunotherapy (Interaction P = 0.0056)")

g173_tn <- copy(g173_df[HR_neg == TRUE])
g173_tn[, Arm_lbl := factor(ifelse(Arm == "durvalumab/olaparib", "TNBC Durva + Ola (n = 21)", "TNBC Control (n = 17)"),
                            levels = c("TNBC Control (n = 17)", "TNBC Durva + Ola (n = 21)"))]
y_d2_top <- max(g173_tn$PPARG, na.rm = TRUE)
d2_span  <- diff(range(g173_tn$PPARG, na.rm = TRUE))
brk_d2 <- data.table(
  Arm_lbl = factor(c("TNBC Control (n = 17)", "TNBC Durva + Ola (n = 21)"), levels = levels(g173_tn$Arm_lbl)),
  y_top   = y_d2_top + 0.08 * d2_span,
  y_tick  = y_d2_top + 0.04 * d2_span,
  y_txt   = y_d2_top + 0.13 * d2_span,
  p_lbl   = c("italic(P) == 0.12", "italic(P) == 0.020")
)

p8d2 <- ggplot(g173_tn, aes(x = Resp, y = PPARG, fill = Resp)) +
  geom_boxplot(width = 0.45, alpha = 0.85, outlier.size = 0.8) +
  geom_jitter(width = 0.13, size = 1.15, alpha = 0.55) +
  geom_segment(data = brk_d2, aes(x = 1, xend = 2, y = y_top, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_segment(data = brk_d2, aes(x = 1, xend = 1, y = y_tick, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_segment(data = brk_d2, aes(x = 2, xend = 2, y = y_tick, yend = y_top), inherit.aes = FALSE, color = "#0F172A", linewidth = 0.55) +
  geom_text(data = brk_d2, aes(x = 1.5, y = y_txt, label = p_lbl), inherit.aes = FALSE, parse = TRUE, size = 3.0, color = "#0F172A") +
  facet_wrap(~ Arm_lbl) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  scale_fill_manual(values = c("RD" = pal_npg[["null"]], "pCR" = pal_npg[["blue"]])) +
  labs(title = "D2. I-SPY2 PPARG in TNBC (Interaction P = 0.049)", x = NULL, y = "PPARG Expression") +
  theme_journal(9.0) + theme(legend.position = "none")

sub_d <- fread(file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
d_dedup <- sub_d[order(gene_symbol, all_spearman_pval)][, .SD[1], by = .(gene_symbol, DRUG_NAME)]
d_top <- d_dedup[all_spearman_pval < 0.05][order(gene_symbol, all_spearman_pval)][, .SD[1:min(9, .N)], by = gene_symbol]
d_top[, drug_clean := gsub("_", "-", DRUG_NAME)]
d_top[, dir := fifelse(subtype_adjusted_pval < 0.05,
               fifelse(all_spearman_rho < 0, "Sensitive (Subtype-Adj. P < 0.05)", "Resistant (Subtype-Adj. P < 0.05)"),
               fifelse(all_spearman_rho < 0, "Sensitive (Unadj. P < 0.05)", "Resistant (Unadj. P < 0.05)"))]
d_top[, drug_ord := paste0(gene_symbol, "___", drug_clean)]
d_top <- d_top[order(gene_symbol, -all_spearman_pval)]; d_top[, drug_ord := factor(drug_ord, levels = drug_ord)]

p8e <- ggplot(d_top, aes(x = all_spearman_rho, y = drug_ord, color = dir, size = -log10(all_spearman_pval))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_segment(aes(x = 0, xend = all_spearman_rho, y = drug_ord, yend = drug_ord), linewidth = 0.65, show.legend = FALSE) +
  geom_point() + facet_wrap(~ gene_symbol, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_x_continuous(breaks = c(-0.50, -0.25, 0, 0.25, 0.50)) +
  coord_cartesian(xlim = c(-0.56, 0.58)) +
  scale_color_manual(values = c(
    "Sensitive (Subtype-Adj. P < 0.05)" = pal_npg[["green"]], "Sensitive (Unadj. P < 0.05)" = "#6EE7B7",
    "Resistant (Subtype-Adj. P < 0.05)" = pal_npg[["red"]],   "Resistant (Unadj. P < 0.05)" = "#FCA5A5"
  )) +
  labs(title = "E. GDSC2 Drug Sensitivity (n = 47 Cell Lines)",
       x = expression(bold(Spearman~rho~"(Expression vs. ln IC"[50]*")")), y = NULL, color = NULL,
       size = expression(bold(-log[10](italic(P))))) +
  theme_journal(8.8) + theme(legend.position = "bottom")

gb_mat <- readRDS(file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))
gb_mat[, Subtype := factor(Subtype_Group, levels = c("ER_pos", "HER2_pos", "TNBC"), labels = c("ER+", "HER2+", "TNBC"))]

p8f1 <- ggplot(gb_mat[DRUG_NAME == "Dasatinib"], aes(x = CXCL8, y = LN_IC50, color = Subtype)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75) +
  geom_smooth(data = gb_mat[DRUG_NAME == "Dasatinib" & Subtype_Group == "TNBC"], method = "lm", se = FALSE, color = pal_npg[["red"]], linewidth = 0.85) +
  scale_color_manual(values = c("ER+" = pal_npg[["blue"]], "HER2+" = pal_npg[["amber"]], "TNBC" = pal_npg[["red"]])) +
  labs(title = "F1. CXCL8 vs. Dasatinib Sensitivity",
       x = expression(bold(CXCL8~(log[2]~TPM + 1))), y = expression(bold(Dasatinib~ln~IC[50])), color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

p8f2 <- ggplot(gb_mat[DRUG_ID == 1200], aes(x = PPARG, y = LN_IC50, color = Subtype)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75) +
  geom_smooth(data = gb_mat[DRUG_ID == 1200 & Subtype_Group == "ER_pos"], method = "lm", se = FALSE, color = pal_npg[["blue"]], linewidth = 0.85) +
  scale_color_manual(values = c("ER+" = pal_npg[["blue"]], "HER2+" = pal_npg[["amber"]], "TNBC" = pal_npg[["red"]])) +
  labs(title = "F2. PPARG vs. Fulvestrant Resistance",
       x = expression(bold(PPARG~(log[2]~TPM + 1))), y = expression(bold(Fulvestrant~ln~IC[50])), color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

fig8_r1 <- (p8a | p8b | km_8c) + plot_layout(widths = c(1.18, 0.94, 0.94))
fig8_r2 <- patchwork::wrap_elements(full = (p8d1 | p8d2))
fig8_r3 <- (p8e | p8f1 | p8f2) + plot_layout(widths = c(1.15, 0.92, 0.92))
fig8_all <- fig8_r1 / fig8_r2 / fig8_r3 + plot_layout(heights = c(1.15, 0.95, 1.05))
save_pub(fig8_all, "Fig8_Therapy_Response_ISPY2_and_GDSC2", w = 18.2, h = 16.5)
rm(fig8_all, fig8_r1, fig8_r2, fig8_r3, p8a, p8b, km_8c, p8d1, p8d2, p8e, p8f1, p8f2, gb_mat); gc()

# ==============================================================================
# FIG S1: 上游因果菌群与代谢物 MR 多方法散点与 LOO 敏感性 (去副标题 + 开放主题)
# ==============================================================================
single_dt <- fread(file.path(res_dir, "planA_mr_single_snp_wald.csv"))
sens_dt   <- fread(file.path(res_dir, "planA_mr_sensitivity_summary.csv"))
loo_dt    <- fread(file.path(res_dir, "planA_mr_leaveoneout.csv"))

target_s1 <- list(
  list(eid = "ebi-a-GCST90016960", oid = "ieu-a-1127", lbl = "Adlercreutzia \u2192 ER+"),
  list(eid = "ebi-a-GCST90016987", oid = "ieu-a-1128", lbl = "Desulfovibrio \u2192 ER\u2212"),
  list(eid = "ebi-a-GCST90017039", oid = "ieu-a-1126", lbl = "Parabacteroides \u2192 Overall"),
  list(eid = "ebi-a-GCST90016951", oid = "ieu-a-1127", lbl = "Ruminococcaceae \u2192 ER+"),
  list(eid = "ebi-a-GCST90016989", oid = "ieu-a-1128", lbl = "Dorea \u2192 ER\u2212"),
  list(eid = "met-d-Acetate",      oid = "ieu-a-1126", lbl = "Acetate \u2192 Overall")
)

sc_plots  <- list()
loo_plots <- list()
for (tk in target_s1) {
  sd <- single_dt[id.exposure == tk$eid & id.outcome == tk$oid]
  flip <- sd$beta.exposure < 0
  sd[, bx := ifelse(flip, -beta.exposure, beta.exposure)]
  sd[, by := ifelse(flip, -beta.outcome, beta.outcome)]
  ms <- sens_dt[id.exposure == tk$eid & id.outcome == tk$oid &
                  method %in% c("Inverse variance weighted", "Weighted median", "MR Egger", "Weighted mode")]
  p_sc <- ggplot(sd, aes(x = bx, y = by)) +
    geom_hline(yintercept = 0, color = "#CBD5E1", linetype = "dotted") +
    geom_errorbar(aes(ymin = by - 1.96 * se.outcome, ymax = by + 1.96 * se.outcome), color = "#94A3B8", width = 0, linewidth = 0.4) +
    geom_errorbar(aes(xmin = bx - 1.96 * se.exposure, xmax = bx + 1.96 * se.exposure), orientation = "y", color = "#94A3B8", width = 0, linewidth = 0.4) +
    geom_abline(data = ms, aes(slope = b, intercept = ifelse(method == "MR Egger" & !is.na(egger_intercept), egger_intercept, 0), color = method), linewidth = 0.85) +
    geom_point(size = 2.0, color = "#0F172A") +
    scale_color_manual(values = c("Inverse variance weighted" = pal_npg[["red"]], "Weighted median" = pal_npg[["green"]],
                                  "MR Egger" = pal_npg[["blue"]], "Weighted mode" = pal_npg[["amber"]])) +
    labs(title = tk$lbl, x = "SNP Effect on Exposure", y = "SNP Effect on BRCA") +
    theme_journal(9.0) + theme(legend.position = "none")
  sc_plots[[tk$lbl]] <- p_sc

  ld <- loo_dt[id.exposure == tk$eid & id.outcome == tk$oid]
  ld[, is_all := (SNP == "All")]
  ld <- ld[order(is_all, -b)]
  ld[, snp_ord := factor(SNP, levels = unique(rev(SNP)))]
  p_lo <- ggplot(ld, aes(x = b, y = snp_ord, color = is_all)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
    geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), orientation = "y", width = 0.25, linewidth = 0.55) +
    geom_point(size = 1.95) +
    scale_color_manual(values = c("FALSE" = "#475569", "TRUE" = pal_npg[["red"]]), guide = "none") +
    labs(title = paste("LOO:", tk$lbl), x = expression(bold("Leave-One-Out IVW"~beta)), y = NULL) +
    theme_journal(9.0)
  loo_plots[[tk$lbl]] <- p_lo
}

fig_s1 <- wrap_plots(c(sc_plots, loo_plots), ncol = 3) +
  plot_annotation(
    title = "Supplementary Figure S1. Upstream MR Multi-Method Sensitivity & Leave-One-Out Analysis",
    theme = theme(plot.title = element_text(face="bold", size=12, family="sans", color="#0F172A"))
  )
save_pub(fig_s1, "FigS1_Upstream_MR_Sensitivity", w = 15.5, h = 14.5)
rm(fig_s1, sc_plots, loo_plots); gc()

# ==============================================================================
# FIG S2: 异质性分解 (洁净Y轴命名 + 正规 I²) + 跨队列欧洲菌群 Meta (无压线)
# ==============================================================================
het_map <- c(
  "All_5_Cohorts_Overall_OS"        = "All 5 Cohorts (OS)",
  "LOO_Excluding_METABRIC"          = "Excl. METABRIC",
  "LOO_Excluding_GSE42568"          = "Excl. GSE42568",
  "LOO_Excluding_TCGA_BRCA"         = "Excl. TCGA-BRCA",
  "LOO_Excluding_GSE20685"          = "Excl. GSE20685",
  "LOO_Excluding_GSE7390"           = "Excl. GSE7390",
  "Microarray_4Cohorts_Overall_OS"  = "Microarray Cohorts Only (n = 4)",
  "GEO_3Cohorts_Affy_OS"            = "Homogeneous Affymetrix GEO (n = 3)",
  "GEO_Large_GSE20685_GSE7390_OS"   = "Large Affymetrix GEO (GSE20685 + GSE7390)",
  "ER_neg_Stratum_TCGA_METABRIC_OS" = "TCGA + METABRIC ER\u2212 Stratum"
)

het_dt <- fread(file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))
het_dt[, lbl_clean := ifelse(analysis %in% names(het_map), unname(het_map[analysis]), gsub("_", " ", analysis))]
het_dt[, lbl := factor(lbl_clean, levels = rev(lbl_clean))]
ps2a <- ggplot(het_dt, aes(x = RE_HR, y = lbl, color = I2_pct == 0)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = RE_LCI, xmax = RE_UCI), orientation = "y", width = 0.26, linewidth = 0.72) +
  geom_point(size = 2.75) +
  geom_text(aes(x = RE_UCI + 0.015, label = sprintf("RE HR = %.3f, I\u00b2 = %.1f%%, P = %.2g", RE_HR, I2_pct, RE_pval)),
            hjust = 0, size = 2.55, color = "#0F172A", show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = pal_npg[["red"]], "FALSE" = pal_npg[["blue"]]),
                     labels = c("TRUE" = "I\u00b2 = 0% (Homogeneous)", "FALSE" = "I\u00b2 > 0% (Heterogeneous)")) +
  coord_cartesian(xlim = c(0.88, 1.72)) +
  labs(title = "A. CXCL8 OS Meta Heterogeneity Audit", x = "Random-Effects HR (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.2) + theme(legend.position = "bottom")

cc_meta  <- fread(file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))
cc_rep   <- fread(file.path(res_dir, "module3_dmp_fr02_core6_replication.csv"))
mib_core <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))[gwas_name != "Marvinbryantia"]

core_pairs <- rbind(
  data.table(pair="Dorea \u2192 ER\u2212", coh="MiBioGen (Discovery)", OR=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR][1], lci=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Dorea \u2192 ER\u2212", coh="DMP (Replication)", OR=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Dorea \u2192 ER\u2212", coh="Pooled FE Meta (k = 2)", OR=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_OR][1], lci=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_LCI][1], uci=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_UCI][1], p=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_pval][1]),
  data.table(pair="Desulfovibrio \u2192 ER\u2212", coh="MiBioGen (Discovery)", OR=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR][1], lci=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio \u2192 ER\u2212", coh="DMP (D. piger)", OR=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio \u2192 ER\u2212", coh="FR02 (D. piger)", OR=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio \u2192 ER\u2212", coh="Pooled FE Meta (k = 3)", OR=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_OR][1], lci=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_LCI][1], uci=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_UCI][1], p=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_pval][1]),
  data.table(pair="Adlercreutzia \u2192 ER+", coh="MiBioGen (Discovery)", OR=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR][1], lci=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_LCI95][1], uci=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_UCI95][1], p=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", pval][1]),
  data.table(pair="Adlercreutzia \u2192 ER+", coh="DMP (Replication)", OR=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR][1], lci=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_LCI95][1], uci=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_UCI95][1], p=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", pval][1]),
  data.table(pair="Adlercreutzia \u2192 ER+", coh="Pooled FE Meta (k = 2)", OR=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_OR][1], lci=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_LCI][1], uci=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_UCI][1], p=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_pval][1])
)
core_pairs[, is_meta := grepl("Meta", coh)]
core_pairs[, row_ord := factor(paste(pair, coh, sep=" | "), levels = rev(paste(pair, coh, sep=" | ")))]

ps2b <- ggplot(core_pairs, aes(x = OR, y = row_ord, color = is_meta)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = lci, xmax = uci), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(aes(shape = is_meta), size = 2.75) +
  geom_text(aes(x = 1.13, label = sprintf("OR = %.3f, P = %.3g", OR, p)), hjust = 0, size = 2.55, color = "#0F172A", show.legend = FALSE) +
  facet_wrap(~ pair, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
  scale_color_manual(values = c("FALSE" = pal_npg[["blue"]], "TRUE" = pal_npg[["red"]]), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  coord_cartesian(xlim = c(0.66, 1.34)) +
  labs(title = "B. Cross-Cohort Gut Microbiome Replication", x = "Odds Ratio (95% CI)", y = NULL) +
  theme_journal(9.2)

fig_s2 <- (ps2a | ps2b) + plot_layout(widths = c(1.05, 1.0))
save_pub(fig_s2, "FigS2_HeterogeneityLOO_and_MicrobiomeMeta", w = 17.0, h = 8.2)
rm(fig_s2, ps2a, ps2b); gc()

# ==============================================================================
# FIG S3: 外部双单细胞验证 + HPA 参考单细胞 + GSE9195 DMFS KM
# ==============================================================================
em_df <- copy(sc_cache[["BRCA_EMTAB8107"]])
ct_em_clean <- clean_lineage_str(em_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_em <- setNames(colorRampPalette(c(pal_npg[["blue"]], pal_npg[["red"]], pal_npg[["green"]], pal_npg[["cyan"]],
                                      pal_npg[["coral"]], pal_npg[["slate"]], pal_npg[["amber"]], "#7C3AED", "#EC4899", "#14B8A6", "#64748B"))(length(ct_em_clean)), ct_em_clean)
ps3a1 <- plot_umap_ct(em_df, pal_em, "A1. E-MTAB-8107 Lineages (n = 33,043)")
ps3a2 <- plot_umap_gene(em_df, "PPARG", "A2. E-MTAB-8107: PPARG", hi_col = pal_npg[["blue"]])
ps3a3 <- plot_umap_gene(em_df, "CXCL8", "A3. E-MTAB-8107: CXCL8", hi_col = pal_npg[["red"]])

g14_df <- copy(sc_cache[["BRCA_GSE148673"]])
ct_g14_clean <- clean_lineage_str(g14_df[, .N, by = Celltype_major][order(-N)]$Celltype_major)
pal_g14 <- setNames(colorRampPalette(c(pal_npg[["blue"]], pal_npg[["red"]], pal_npg[["green"]], pal_npg[["cyan"]],
                                       pal_npg[["coral"]], pal_npg[["slate"]], pal_npg[["amber"]], "#7C3AED", "#14B8A6"))(length(ct_g14_clean)), ct_g14_clean)
ps3b1 <- plot_umap_ct(g14_df, pal_g14, "B1. GSE148673 Lineages (n = 10,359)")
ps3b2 <- plot_umap_gene(g14_df, "PPARG", "B2. GSE148673: PPARG", hi_col = pal_npg[["blue"]])
ps3b3 <- plot_umap_gene(g14_df, "IL1B",  "B3. GSE148673: IL1B",  hi_col = pal_npg[["red"]])

sc_hpa <- fread(file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
sc_keep <- c("Neutrophils","Macrophages","Fibroblasts","Adipocytes","Lymphatic endothelial cells","Vascular endothelial cells","Breast secretory cells","T-cells","B-cells")
sc_sub <- sc_hpa[gene_symbol %in% c("CXCL8","PPARG") & cell_type %in% sc_keep]
sc_sub[, cell_clean := tools::toTitleCase(gsub("-", " ", cell_type))]
sc_sub[, log10_ncpm := log10(nCPM + 1)]
sc_ord <- sc_sub[, .(mx = max(log10_ncpm)), by = cell_clean][order(mx)]$cell_clean
sc_sub[, cell_clean := factor(cell_clean, levels = sc_ord)]
ps3c <- ggplot(sc_sub, aes(x = log10_ncpm, y = cell_clean, fill = gene_symbol)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.64) +
  scale_fill_manual(values = c("CXCL8" = pal_npg[["red"]], "PPARG" = pal_npg[["blue"]])) +
  labs(title = "C. HPA Reference Single-Cell Expression", x = expression(bold(log[10](nCPM + 1))), y = NULL, fill = NULL) +
  theme_journal(9.0) + theme(legend.position = "bottom")

g91_km <- copy(th_pc$gse9195)[!is.na(t_dmfs) & t_dmfs > 0 & !is.na(e_dmfs)]
g91_km[, Group := ifelse(PPARG >= median(PPARG), "High", "Low")]
km_s3d <- plot_km_with_risk(
  g91_km, "t_dmfs", "e_dmfs", "Group",
  "D. GSE9195 Tamoxifen ER+ DMFS", "Years", "Distant Metastasis-Free Survival", "PPARG",
  c("Low", "High"), c(pal_npg[["red"]], pal_npg[["blue"]]),
  pval_coord = c(0.4, 0.16), break_by = 3, max_x = 12, base_size = 9.0
)

figs3_r1 <- (ps3a1 | ps3a2 | ps3a3)
figs3_r2 <- (ps3b1 | ps3b2 | ps3b3)
figs3_r3 <- patchwork::wrap_elements(full = (ps3c | km_s3d) + plot_layout(widths = c(1.05, 1.05)))
fig_s3 <- figs3_r1 / figs3_r2 / figs3_r3 + plot_layout(heights = c(0.95, 0.95, 1.15))
save_pub(fig_s3, "FigS3_scRNA_ExternalCohorts_and_GSE9195", w = 17.2, h = 16.2)

cat("\nALL 11 NATURE/CELL/FRONTIERS FIGURES COMPLETE!\n")
