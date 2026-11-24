# 14_generate_draft_figures.R
# 生成模块 1~6 初版高质量科研图表 (保存 PDF + 300dpi PNG 至 03_figures/)

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(survminer)
})

res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
fig_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures")
bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F4F8", color = "#CBD5E1"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.2, color = "#475569")
  )

save_both <- function(p, name, w = 12, h = 7) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  cat("Saved figure:", name, "\n")
}

# ==============================================================================
# Fig 1: 上游因果菌群 (6个有效菌群) & 血浆代谢物 (5个有效代谢物) MR 森林图
# ==============================================================================
mib_sig <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))
mib_sig <- mib_sig[gwas_name != "Marvinbryantia"]
mib_sig[, label := paste0(gwas_rank, " ", gwas_name, " -> ", outcome_label)]
mib_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", fifelse(FDR < 0.10, "FDR < 0.10", "Nominal P < 0.05"))]
mib_sig <- mib_sig[order(pval)]
mib_sig[, label := factor(label, levels = rev(label))]

p1a <- ggplot(mib_sig, aes(x = OR, y = label, color = sig_tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(x = max(OR_UCI95) * 1.04, label = sprintf("OR=%.2f, P=%.2g", OR, pval)),
            hjust = 0, size = 3, color = "#1E293B", show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = "#DC2626", "FDR < 0.10" = "#EA580C", "Nominal P < 0.05" = "#2563EB")) +
  coord_cartesian(xlim = c(0.65, max(mib_sig$OR_UCI95) * 1.38)) +
  labs(title = "A. Causal Gut Microbiota -> BCAC Breast Cancer (Reverse-MR Filtered)",
       subtitle = "MiBioGen IVW MR (n = 6 valid taxa; Marvinbryantia excluded by reverse MR)",
       x = "Odds Ratio (95% CI)", y = NULL, color = "Significance") +
  theme_pub + theme(legend.position = "bottom")

met_sig <- fread(file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv"))
met_sig[, label := paste0(trait_name, " -> ", outcome_label)]
met_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05")]
met_sig <- met_sig[order(pval)]
met_sig[, label := factor(label, levels = rev(label))]

p1b <- ggplot(met_sig, aes(x = OR, y = label, color = sig_tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(x = max(OR_UCI95) * 1.04, label = sprintf("OR=%.2f, P=%.2g", OR, pval)),
            hjust = 0, size = 3, color = "#1E293B", show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = "#DC2626", "Nominal P < 0.05" = "#2563EB")) +
  coord_cartesian(xlim = c(0.45, max(met_sig$OR_UCI95) * 1.38)) +
  labs(title = "B. Causal Plasma Metabolites -> BCAC Breast Cancer",
       subtitle = "47 gutMGene-aligned metabolites (all 5 hits passed reverse MR, P > 0.05)",
       x = "Odds Ratio (95% CI)", y = NULL, color = "Significance") +
  theme_pub + theme(legend.position = "bottom")

save_both(p1a | p1b, "Fig1_Draft_Microbe_Metabolite_MR", w = 15.5, h = 6.5)

# ==============================================================================
# Fig 2: gutMGene 三元调控网络图
# ==============================================================================
edges_dt <- fread(file.path(res_dir, "module3_candidate_gutmgene_edges.csv"))
cis_pos  <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_g    <- unique(cis_pos$gene_symbol)

sub_e <- edges_dt[gene_symbol %in% pos_g & !is.na(metabolite)]
map_microbe <- function(m) {
  if (grepl("Parabacteroides", m, ignore.case = TRUE)) return("Parabacteroides (Prot.)")
  if (grepl("Dorea", m, ignore.case = TRUE)) return("Dorea (ER- Prot.)")
  if (grepl("Ruminococc", m, ignore.case = TRUE)) return("Ruminococcaceae (Risk)")
  if (grepl("Blautia", m, ignore.case = TRUE)) return("Blautia (ER- Risk)")
  if (grepl("MR_Positive_Metabolite", m)) return("Circulating Acetate (Risk)")
  return(NA_character_)
}
sub_e[, m_grp := sapply(microbe, map_microbe)]
mm_links <- unique(sub_e[!is.na(m_grp), .(from = m_grp, to = metabolite)])
mm_links <- unique(rbind(mm_links, data.table(from = "Adlercreutzia (ER+ Prot.)", to = "Equol")))
mg_links <- unique(sub_e[, .(from = metabolite, to = gene_symbol, alteration)])
mg_links <- mg_links[, .(alteration = ifelse(length(unique(alteration)) > 1, "mixed", alteration[1])), by = .(from, to)]

left_nodes  <- sort(unique(mm_links$from))
mid_nodes   <- sort(unique(c(mm_links$to, mg_links$from)))
right_nodes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "HDAC11", "TNF", "IL12A", "CXCR2")

node_pos <- rbindlist(list(
  data.table(node = left_nodes,  col = 1, y = seq(1, 10, length.out = length(left_nodes)),  layer = "Causal Microbe / Blood"),
  data.table(node = mid_nodes,   col = 2, y = seq(1.5, 9.5, length.out = length(mid_nodes)), layer = "gutMGene Metabolite"),
  data.table(node = right_nodes, col = 3, y = seq(1, 10, length.out = length(right_nodes)), layer = "cis-eQTL Host Gene")
))

e1 <- merge(mm_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e1 <- merge(e1, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to")
e1[, edge_type := "produces"]

e2 <- merge(mg_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e2 <- merge(e2, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to")
e2[, edge_type := alteration]

all_seg <- rbindlist(list(e1[, .(x1, y1, x2, y2, edge_type)], e2[, .(x1, y1, x2, y2, edge_type)]))
node_pos[, node_fill := fifelse(col == 1, "#E0F2FE",
                        fifelse(col == 2, "#FEF3C7",
                        fifelse(node %in% c("PPARG", "CXCL8"), "#DCFCE7",
                        fifelse(node == "CXCR2", "#FEE2E2", "#F1F5F9"))))]

p2 <- ggplot() +
  geom_segment(data = all_seg, aes(x = x1 + 0.14, y = y1, xend = x2 - 0.14, yend = y2, color = edge_type),
               linewidth = 0.85, alpha = 0.78) +
  geom_label(data = node_pos, aes(x = col, y = y, label = node, fill = node_fill),
             size = 3.4, fontface = "bold", label.padding = unit(0.28, "lines"), color = "#0F172A") +
  scale_fill_identity() +
  scale_color_manual(
    values = c("produces" = "#64748B", "activation" = "#16A34A", "inhibition" = "#DC2626", "mixed" = "#9333EA"),
    name = "gutMGene Regulation"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("Causal Gut Microbes / Metabolites (MR)", "Microbial Metabolites (gutMGene)", "Host Target Genes (cis-eQTL MR)"),
                     limits = c(0.6, 3.45)) +
  labs(title = "Gut Microbiota - Metabolite - Host Gene Causal Regulatory Network",
       subtitle = "Green host genes: top validated stars (PPARG ER+, CXCL8 ER-); Red: CXCR2 excluded by coloc PP.H3=0.9996",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

save_both(p2, "Fig2_Draft_gutMGene_Causal_Network", w = 13, h = 7)

# ==============================================================================
# Fig 3: 8 个纯 cis-eQTL 阳性基因 MR 森林图 + coloc 共定位后验概率柱状图
# ==============================================================================
cis_pos[, label := paste0(gene_symbol, " -> ", outcome_label)]
cis_pos[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05")]
cis_pos <- cis_pos[order(pval)]
cis_pos[, label := factor(label, levels = rev(label))]

p3a <- ggplot(cis_pos, aes(x = OR, y = label, color = sig_tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(aes(x = max(OR_UCI95) * 1.04, label = sprintf("OR=%.2f, P=%.2g", OR, pval)),
            hjust = 0, size = 2.9, color = "#1E293B", show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = "#DC2626", "Nominal P < 0.05" = "#2563EB")) +
  coord_cartesian(xlim = c(0.45, max(cis_pos$OR_UCI95) * 1.4)) +
  labs(title = "A. Strict cis-eQTLGen -> BCAC MR (8 Genes, 13 Pairs)",
       subtitle = "Filtered from 49 candidate genes (33 with cis-IVs within +-1 Mb)",
       x = "Odds Ratio (95% CI)", y = NULL, color = "Significance") +
  theme_pub + theme(legend.position = "bottom")

coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
coloc_dt[, pair := paste0(gene_symbol, " (", outcome_label, ")")]
coloc_long <- melt(coloc_dt, id.vars = "pair", measure.vars = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
                   variable.name = "Hypothesis", value.name = "PP")
coloc_dt <- coloc_dt[order(PP.H4)]
coloc_long[, pair := factor(pair, levels = coloc_dt$pair)]

p3b <- ggplot(coloc_long, aes(x = PP, y = pair, fill = Hypothesis)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("PP.H0"="#CBD5E1","PP.H1"="#94A3B8","PP.H2"="#64748B","PP.H3"="#F97316","PP.H4"="#16A34A")) +
  labs(title = "B. Bayesian Colocalization (eQTLGen vs BCAC, +-500kb)",
       subtitle = "IL1B PP.H4=0.507 (H4/[H3+H4]=0.895); CXCR2 flagged LD-confounded (H3=0.9996)",
       x = "Posterior Probability", y = NULL) +
  theme_pub + theme(legend.position = "bottom")

save_both(p3a | p3b, "Fig3_Draft_cisEQTL_MR_and_Coloc", w = 15.5, h = 6.8)

# ==============================================================================
# Fig 4: 核心明星基因 PPARG (ER+) 与 CXCL8 (ER-) 多队列差异表达与生存图
# ==============================================================================
tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]

df_pparg_tn <- rbind(
  data.table(Group = "Tumor (n=1082)", Expr = as.numeric(tc$tumor_log2["PPARG", ])),
  data.table(Group = "Normal (n=114)", Expr = as.numeric(tc$normal_log2["PPARG", ]))
)
df_pparg_tn[, Group := factor(Group, levels = c("Normal (n=114)", "Tumor (n=1082)"))]
p4a <- ggplot(df_pparg_tn, aes(x = Group, y = Expr, fill = Group)) +
  geom_violin(alpha = 0.45, trim = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Normal (n=114)" = "#38BDF8", "Tumor (n=1082)" = "#F87171")) +
  labs(title = "A. TCGA-BRCA: PPARG Tumor vs Normal",
       subtitle = "log2FC = -3.23, Wilcoxon P = 7.74e-57", x = NULL, y = "log2(RSEM + 1)") +
  theme_pub + theme(legend.position = "none")

idx_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_km_p <- data.table(
  time = cl_tc$OS_MONTHS[idx_erp],
  status = cl_tc$OS_event[idx_erp],
  expr = as.numeric(tc$tumor_log2["PPARG", idx_erp])
)
df_km_p[, Group := factor(ifelse(expr >= median(expr), "High PPARG", "Low PPARG"), levels = c("Low PPARG", "High PPARG"))]
fit_p <- survfit(Surv(time, status) ~ Group, data = df_km_p)
sp_p <- survminer::ggsurvplot(fit_p, data = df_km_p, pval = TRUE, risk.table = FALSE,
                              palette = c("#EF4444", "#2563EB"),
                              title = "B. TCGA ER+ OS: PPARG (n=783)",
                              xlab = "Time (Months)", ylab = "Overall Survival",
                              ggtheme = theme_pub)
p4b <- sp_p$plot + labs(subtitle = "Cox HR = 0.820, P = 0.0135 | Log-rank P = 0.00084")

mb <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
df_cx_er <- rbind(
  data.table(Cohort = "TCGA-BRCA", ER = ifelse(cl_tc$ER_Xena == "Positive", "ER+", ifelse(cl_tc$ER_Xena == "Negative", "ER-", NA_character_)),
             Expr = as.numeric(tc$tumor_log2["CXCL8", ])),
  data.table(Cohort = "METABRIC", ER = ifelse(mb$clin$ER_clean == "Positive", "ER+", ifelse(mb$clin$ER_clean == "Negative", "ER-", NA_character_)),
             Expr = as.numeric(mb$expr_target["CXCL8", ]))
)
df_cx_er <- df_cx_er[!is.na(ER)]
p4c <- ggplot(df_cx_er, aes(x = ER, y = Expr, fill = ER)) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c("ER+" = "#60A5FA", "ER-" = "#FB7185")) +
  labs(title = "C. CXCL8 Expression in ER- vs ER+ Subtypes",
       subtitle = "TCGA P = 9.99e-28 (log2FC=-1.74); METABRIC P = 7.29e-53 (log2FC=-1.20)",
       x = NULL, y = "Normalized log2 Expression") +
  theme_pub + theme(legend.position = "none")

g20 <- readRDS(file.path(geo_dir, "GSE20685_parsed.rds"))
pd20 <- as.data.table(g20$pheno)
df_km_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_km_c[, Group := factor(ifelse(expr >= median(expr), "High CXCL8", "Low CXCL8"), levels = c("Low CXCL8", "High CXCL8"))]
fit_c <- survfit(Surv(time, status) ~ Group, data = df_km_c)
sp_c <- survminer::ggsurvplot(fit_c, data = df_km_c, pval = TRUE, risk.table = FALSE,
                              palette = c("#2563EB", "#EF4444"),
                              title = "D. GSE20685 OS: CXCL8 (n=327)",
                              xlab = "Time (Years)", ylab = "Overall Survival",
                              ggtheme = theme_pub)
p4d <- sp_c$plot + labs(subtitle = "Cox HR = 1.251, P = 0.00057 | Log-rank P = 0.0160")

# 4E: 修复 stratum 标签以区分 Overall 与 ER_pos 多因素 Cox
surv_all <- fread(file.path(res_dir, "module4_survival_analysis.csv"))
mv_cox   <- fread(file.path(res_dir, "module4_multivariable_cox.csv"))
hr_sub <- rbind(
  surv_all[gene_symbol == "PPARG" & ((cohort == "TCGA_BRCA" & stratum %in% c("Overall","ER_pos") & endpoint == "OS") |
                                     (cohort == "METABRIC" & stratum %in% c("Overall","ER_pos") & endpoint %in% c("OS","RFS")) |
                                     (cohort %in% c("GSE42568","GSE20685","GSE7390") & endpoint == "OS")),
           .(gene_symbol, label = paste(cohort, stratum, endpoint), HR, HR_LCI95, HR_UCI95, pval = cox_pval)],
  surv_all[gene_symbol == "CXCL8" & ((cohort == "METABRIC" & stratum == "Overall" & endpoint %in% c("OS","RFS")) |
                                     (cohort %in% c("GSE42568","GSE20685","GSE7390") & endpoint %in% c("OS","MFS","RFS"))),
           .(gene_symbol, label = paste(cohort, stratum, endpoint), HR, HR_LCI95, HR_UCI95, pval = cox_pval)],
  mv_cox[gene_symbol == "CXCL8" & cohort == "METABRIC" & stratum == "Overall",
         .(gene_symbol, label = "METABRIC Overall OS_Multivar(Age+Stage)", HR, HR_LCI95, HR_UCI95, pval)]
)
hr_sub[, label := paste0(gene_symbol, " | ", label)]
hr_sub[, sig := ifelse(pval < 0.05, "P < 0.05", "P >= 0.05")]
hr_sub <- hr_sub[order(gene_symbol, pval)]
hr_sub[, label := factor(label, levels = rev(unique(label)))]

p4e <- ggplot(hr_sub, aes(x = HR, y = label, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.25, linewidth = 0.65) +
  geom_point(size = 2.6) +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("P < 0.05" = "#DC2626", "P >= 0.05" = "#64748B")) +
  labs(title = "E. Multi-Cohort Prognostic Consistency Forest Plot",
       subtitle = "PPARG protects ER+ BRCA (HR<1); CXCL8 drives poor prognosis across independent cohorts (HR>1)",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_pub + theme(legend.position = "bottom")

rm(tc, mb, g20); gc()
fig4 <- (p4a | p4b) / (p4c | p4d)
save_both(fig4, "Fig4_Draft_PPARG_CXCL8_Expr_KM", w = 13.5, h = 10)
save_both(p4e, "Fig4E_Draft_MultiCohort_Survival_Forest", w = 10, h = 7.5)

# ==============================================================================
# Fig 5: 模块 5 免疫微环境相关性 & Hallmark GSEA (分面内独立排序)
# ==============================================================================
cor5 <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
tme_top <- cor5[((gene_symbol == "CXCL8" & stratum == "ER_neg") | (gene_symbol == "PPARG" & stratum == "ER_pos")) &
                  feature_type == "TME_Cell_Score" & FDR < 0.05]
tme_top <- tme_top[order(gene_symbol, -abs(rho))][, .SD[1:min(10, .N)], by = gene_symbol]
tme_top[, facet_lbl := ifelse(gene_symbol == "CXCL8", "CXCL8 in ER- BRCA (n=237)", "PPARG in ER+ BRCA (n=783)")]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feature)]
tme_top <- tme_top[order(facet_lbl, rho)]
tme_top[, feat_ord := factor(feat_ord, levels = feat_ord)]

p5a <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_segment(aes(x = 0, xend = rho, y = feat_ord, yend = feat_ord), linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("Positive" = "#DC2626", "Negative" = "#2563EB")) +
  labs(title = "A. Tumor Microenvironment Correlations (FDR < 0.05)",
       x = "Spearman rho", y = NULL, color = NULL) +
  theme_pub + theme(legend.position = "bottom")

gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene == "CXCL8" & stratum == "ER_neg") | (target_gene == "PPARG" & stratum == "ER_pos")) & p.adjust < 0.05]
gs_top <- gs_top[order(target_gene, -abs(NES))][, .SD[1:min(10, .N)], by = target_gene]
gs_top[, clean_id := gsub("^HALLMARK_", "", ID)]
gs_top[, facet_lbl := ifelse(target_gene == "CXCL8", "CXCL8 High vs Low (ER- BRCA)", "PPARG High vs Low (ER+ BRCA)")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]
gs_top[, id_ord := factor(id_ord, levels = id_ord)]

p5b <- ggplot(gs_top, aes(x = NES, y = id_ord, fill = -log10(p.adjust))) +
  geom_col(width = 0.72) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_fill_gradient(low = "#93C5FD", high = "#1D4ED8") +
  labs(title = "B. Single-Gene Hallmark GSEA (FDR < 0.05)",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = "-log10(FDR)") +
  theme_pub + theme(legend.position = "bottom")

save_both(p5a / p5b, "Fig5_Draft_TME_and_Hallmark_GSEA", w = 13.5, h = 10)

# ==============================================================================
# Fig 6: 模块 6 GDSC2 乳腺癌细胞系实测药敏关联 (按药名去重保留最强记录)
# ==============================================================================
drug6 <- fread(file.path(res_dir, "module6_gdsc2_all_drugs_correlation.csv"))
d_sub <- drug6[gene_symbol %in% c("PPARG", "CXCL8", "IL1B", "FFAR2") & pval_IC50 < 0.05]
d_sub <- d_sub[order(gene_symbol, pval_IC50)][, .SD[1], by = .(gene_symbol, DRUG_NAME)]
d_sub <- d_sub[order(gene_symbol, pval_IC50)][, .SD[1:min(8, .N)], by = gene_symbol]
d_sub[, dir := ifelse(spearman_rho_IC50 < 0, "More Sensitive (rho < 0)", "More Resistant (rho > 0)")]
d_sub[, drug_ord := paste0(gene_symbol, "___", DRUG_NAME)]
d_sub <- d_sub[order(gene_symbol, -pval_IC50)]
d_sub[, drug_ord := factor(drug_ord, levels = drug_ord)]

p6 <- ggplot(d_sub, aes(x = spearman_rho_IC50, y = drug_ord, color = dir, size = -log10(pval_IC50))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_point() +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("More Sensitive (rho < 0)" = "#16A34A", "More Resistant (rho > 0)" = "#DC2626")) +
  labs(title = "GDSC2 Experimental Drug Sensitivity in 47 Breast Cancer Cell Lines",
       subtitle = "Spearman correlation between DepMap24Q4 expression and GDSC2 LN_IC50 (Top P < 0.05 hits per gene)",
       x = "Spearman rho (Gene Expr vs LN_IC50)", y = NULL, color = "Response", size = "-log10(P)") +
  theme_pub

save_both(p6, "Fig6_Draft_GDSC2_Pharmacogenomics", w = 13.5, h = 8.5)
cat("Refined all draft figures cleanly!\n")
