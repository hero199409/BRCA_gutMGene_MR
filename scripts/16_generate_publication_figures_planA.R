# 16_generate_publication_figures_planA.R
# 方案 A（CXCL8 ER- 风险轴 + PPARG ER+ 保护轴）投稿级高清组图生成脚本
# 输出目录: 03_figures/Publication_PlanA/ (300 dpi PNG + PDF)

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(survminer)
  library(jpeg)
  library(grid)
  library(scales)
})

cache_dir <- "00_data/mr_cache"
bulk_dir  <- "00_data/bulk_cohorts"
geo_dir   <- "00_data/geo_cache"
res_dir   <- "02_results"
pub_dir   <- "03_figures/Publication_PlanA"
ihc_dir   <- "03_figures/HPA_IHC_Images"
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

# 统一 NPG 科研期刊配色与主题
pal_npg <- c(
  red   = "#E64B35",
  blue  = "#3C5488",
  cyan  = "#4DBBD5",
  green = "#00A087",
  coral = "#F39B7F",
  slate = "#8491B4",
  amber = "#D97706"
)

theme_journal <- function(base_size = 10.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#F1F5F9", linewidth = 0.4),
      panel.border = element_rect(color = "#334155", linewidth = 0.7),
      strip.background = element_rect(fill = "#F8FAFC", color = "#CBD5E1"),
      strip.text = element_text(face = "bold", size = base_size - 0.5, color = "#0F172A"),
      plot.title = element_text(face = "bold", size = base_size + 1.5, color = "#0F172A"),
      plot.subtitle = element_text(size = base_size - 1.2, color = "#475569"),
      axis.title = element_text(face = "bold", size = base_size - 0.5, color = "#1E293B"),
      axis.text = element_text(color = "#334155", size = base_size - 1.2),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1.2),
      plot.tag = element_text(face = "bold", size = base_size + 3, color = "#0F172A")
    )
}

save_pub <- function(p, name, w = 14, h = 7.5) {
  ggsave(file.path(pub_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(pub_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  cat("Exported:", name, "\n")
}

# 辅助函数：构建森林图通用模块
plot_forest_panel <- function(dt, title_txt, sub_txt, xlim_vals, color_map, text_x_mult = 1.04) {
  ggplot(dt, aes(x = OR, y = label, color = sig_tag)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.6) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.75) +
    geom_point(aes(shape = domain_flag), size = 3.0) +
    geom_text(aes(x = max(OR_UCI95, na.rm = TRUE) * text_x_mult, label = anno_txt),
              hjust = 0, size = 2.9, color = "#1E293B", show.legend = FALSE) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = c("Primary" = 16, "Warning (2 SNPs)" = 17, "LD-Confounded (Excluded)" = 4), guide = "none") +
    coord_cartesian(xlim = xlim_vals) +
    labs(title = title_txt, subtitle = sub_txt, x = "Odds Ratio (95% CI)", y = NULL, color = "Evidence Tier") +
    theme_journal() +
    theme(legend.position = "bottom")
}

# ==============================================================================
# Fig 1: 上游菌群与代谢物两样本 MR (含反向 MR 剔除标识与 Meta 注释)
# ==============================================================================
mib_sig <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))
mib_sig <- mib_sig[gwas_name != "Marvinbryantia"]
mib_sig[, label := sprintf("%s %s -> %s (nSNP=%d)", gwas_rank, gwas_name, gsub("BRCA_", "", outcome_label), nsnp)]
mib_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", fifelse(FDR < 0.10, "FDR < 0.10", "Nominal P < 0.05"))]
mib_sig[, domain_flag := fifelse(nsnp <= 2, "Warning (2 SNPs)", "Primary")]
mib_sig[, anno_txt := sprintf("OR=%.2f [%.2f-%.2f], P=%.2g", OR, OR_LCI95, OR_UCI95, pval)]
mib_sig <- mib_sig[order(pval)]
mib_sig[, label := factor(label, levels = rev(label))]

p1a <- plot_forest_panel(
  mib_sig,
  "A. Causal Gut Microbiota -> Breast Cancer (BCAC)",
  "6 reverse-MR-validated taxa from 64 MiBioGen exposures (Marvinbryantia excluded, rev P=0.022)",
  c(0.65, max(mib_sig$OR_UCI95) * 1.48),
  c("FDR < 0.05" = pal_npg[["red"]], "FDR < 0.10" = pal_npg[["amber"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

met_sig <- fread(file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv"))
# 仅高亮通过 script 07 反向 MR 的 5 大优先代谢物；保留第 6 个关键子串匹配结果作灰色附注或聚焦优先 5 个
met_sig <- met_sig[id.exposure %in% c("met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-446")]
met_sig[, trait_short := fifelse(grepl("hydroxyphen", trait_name), "3-(4-hydroxyphenyl)lactate",
                          fifelse(grepl("phenylpropionate", trait_name), "3-phenylpropionate",
                          fifelse(grepl("linoleic", trait_name), "18:2 Linoleic acid", trait_name)))]
met_sig[, label := sprintf("%s -> %s (nSNP=%d)", trait_short, gsub("BRCA_", "", outcome_label), nsnp)]
met_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05")]
met_sig[, domain_flag := "Primary"]
met_sig[, anno_txt := sprintf("OR=%.2f [%.2f-%.2f], P=%.2g", OR, OR_LCI95, OR_UCI95, pval)]
met_sig <- met_sig[order(pval)]
met_sig[, label := factor(label, levels = rev(label))]

p1b <- plot_forest_panel(
  met_sig,
  "B. Causal Plasma Metabolites -> Breast Cancer (BCAC)",
  "5 prioritized gutMGene metabolites (all 15 reverse-MR pairs P > 0.05; Acetate FDR < 0.05)",
  c(0.45, max(met_sig$OR_UCI95) * 1.48),
  c("FDR < 0.05" = pal_npg[["red"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

save_pub(p1a | p1b, "Fig1_Pub_Upstream_Causal_MR", w = 16.5, h = 6.8)

# ==============================================================================
# Fig 2: 方案 A 双主轴多组学调控网络图（精确区分激活/抑制与排除对照）
# ==============================================================================
edges_dt <- fread(file.path(res_dir, "module3_candidate_gutmgene_edges.csv"))
cis_pos  <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_g    <- unique(cis_pos$gene_symbol)

sub_e <- edges_dt[gene_symbol %in% pos_g & !is.na(metabolite)]
map_microbe <- function(m) {
  if (grepl("Parabacteroides", m, ignore.case = TRUE)) return("Parabacteroides (Prot.)")
  if (grepl("Dorea", m, ignore.case = TRUE)) return("Dorea (ER- Prot.)")
  if (grepl("Ruminococc", m, ignore.case = TRUE)) return("Ruminococcaceae (Risk)")
  if (grepl("Blautia", m, ignore.case = TRUE)) return("Blautia (ER- Risk, 2 SNPs)")
  if (grepl("MR_Positive_Metabolite", m)) return("Blood Acetate (Risk, FDR<0.05)")
  return(NA_character_)
}
sub_e[, m_grp := sapply(microbe, map_microbe)]
mm_links <- unique(sub_e[!is.na(m_grp), .(from = m_grp, to = metabolite)])
mm_links <- unique(rbind(
  mm_links,
  data.table(from = "Adlercreutzia (ER+ Prot., FDR<0.05)", to = "Equol"),
  data.table(from = "Blood Acetate (Risk, FDR<0.05)", to = "Acetate"),
  data.table(from = "Blood Linoleic Acid (MR, P=0.031)", to = "10-keto-12Z-octadecenoic acid")
))
# Keep key metabolites clearly tied to Plan A
mg_links <- unique(sub_e[, .(from = metabolite, to = gene_symbol, alteration)])
mg_links <- mg_links[, .(alteration = ifelse(length(unique(alteration)) > 1, "mixed", alteration[1])), by = .(from, to)]
# Add verified gutMGene PMID:25749343 activation link for PPARG
mg_links <- unique(rbind(
  mg_links,
  data.table(from = "10-keto-12Z-octadecenoic acid", to = "PPARG", alteration = "activation")
))

# Clean ordering along subtype / role hierarchy
left_nodes  <- c("Parabacteroides (Prot.)", "Dorea (ER- Prot.)", "Blautia (ER- Risk, 2 SNPs)",
                 "Blood Acetate (Risk, FDR<0.05)", "Ruminococcaceae (Risk)",
                 "Blood Linoleic Acid (MR, P=0.031)", "Adlercreutzia (ER+ Prot., FDR<0.05)")
mid_nodes   <- c("Butyrate", "Acetate", "Succinate", "3-Indolepropionic acid",
                 "Trimethylamine oxide", "10-keto-12Z-octadecenoic acid", "Equol")
right_nodes <- c("CXCL8 [Plan A ER- Axis]", "IL1B [Coloc PP.H4=0.507]", "FFAR2 [GTEx Breast]",
                 "TNF [Secondary]", "HDAC11 [Secondary]", "IL12A [Secondary]",
                 "PPARG [Plan A ER+ Axis]", "CXCR2 [Excluded: PP.H3=1.0]")

rename_g <- c(
  "CXCL8"  = "CXCL8 [Plan A ER- Axis]",
  "IL1B"   = "IL1B [Coloc PP.H4=0.507]",
  "FFAR2"  = "FFAR2 [GTEx Breast]",
  "TNF"    = "TNF [Secondary]",
  "HDAC11" = "HDAC11 [Secondary]",
  "IL12A"  = "IL12A [Secondary]",
  "PPARG"  = "PPARG [Plan A ER+ Axis]",
  "CXCR2"  = "CXCR2 [Excluded: PP.H3=1.0]"
)
mg_links[, to := rename_g[to]]

# Filter only valid nodes in our curated network
mm_links <- mm_links[from %in% left_nodes & to %in% mid_nodes]
mg_links <- mg_links[from %in% mid_nodes & to %in% right_nodes]

node_pos <- rbindlist(list(
  data.table(node = left_nodes,  col = 1, y = seq(9.5, 1.5, length.out = length(left_nodes)),  layer = "1. Causal Microbe / Metabolite Arm"),
  data.table(node = mid_nodes,   col = 2, y = seq(9.5, 1.5, length.out = length(mid_nodes)),   layer = "2. gutMGene Bioactive Mediators"),
  data.table(node = right_nodes, col = 3, y = seq(9.5, 1.5, length.out = length(right_nodes)), layer = "3. Strict cis-eQTL Host Targets")
))

e1 <- merge(mm_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e1 <- merge(e1, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to")
e1[, edge_type := "produces"]

e2 <- merge(mg_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e2 <- merge(e2, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to")
e2[, edge_type := alteration]

all_seg <- rbindlist(list(e1[, .(x1, y1, x2, y2, edge_type)], e2[, .(x1, y1, x2, y2, edge_type)]))
node_pos[, node_fill := fifelse(grepl("CXCL8", node), "#FFE4E6",
                        fifelse(grepl("PPARG", node), "#DBEAFE",
                        fifelse(grepl("Excluded", node), "#F1F5F9",
                        fifelse(col == 1, "#E0F2FE",
                        fifelse(col == 2, "#FEF3C7", "#ECFDF5")))))]
node_pos[, node_col  := fifelse(grepl("CXCL8", node), "#BE123C",
                        fifelse(grepl("PPARG", node), "#1D4ED8",
                        fifelse(grepl("Excluded", node), "#94A3B8", "#0F172A")))]

p2 <- ggplot() +
  geom_segment(data = all_seg, aes(x = x1 + 0.16, y = y1, xend = x2 - 0.16, yend = y2, color = edge_type, linetype = edge_type),
               linewidth = 0.85, alpha = 0.82) +
  geom_label(data = node_pos, aes(x = col, y = y, label = node, fill = node_fill, color = node_col),
             size = 3.25, fontface = "bold", label.padding = unit(0.32, "lines"), linewidth = 0.6) +
  scale_fill_identity() +
  scale_color_manual(
    values = c("produces" = "#64748B", "activation" = "#059669", "inhibition" = "#DC2626", "mixed" = "#7C3AED",
               "#BE123C" = "#BE123C", "#1D4ED8" = "#1D4ED8", "#94A3B8" = "#94A3B8", "#0F172A" = "#0F172A"),
    breaks = c("produces", "activation", "inhibition", "mixed"),
    name = "gutMGene Edge"
  ) +
  scale_linetype_manual(
    values = c("produces" = "solid", "activation" = "solid", "inhibition" = "dashed", "mixed" = "dotted"),
    breaks = c("produces", "activation", "inhibition", "mixed"),
    name = "gutMGene Edge"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("Causal Microbe / Blood Arm (MR)", "gutMGene Bioactive Metabolites", "Strict cis-eQTL Host Genes (MR)"),
                     limits = c(0.48, 3.55)) +
  ylim(0.8, 10.2) +
  labs(title = "Plan A Multi-Omics Causal Architecture: Subtype-Divergent Microbiota-Metabolite-Host Axes",
       subtitle = "ER- Axis (Rose): Parabacteroides/Dorea -> Butyrate/Acetate -| CXCL8; ER+ Axis (Blue): Linoleic acid -> 10-keto-12Z-octadecenoic acid (PMID:25749343) -> PPARG",
       x = NULL, y = NULL) +
  theme_journal() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

save_pub(p2, "Fig2_Pub_PlanA_DualAxis_Network", w = 14.2, h = 7.4)

# ==============================================================================
# Fig 3: 严格顺式 MR 森林图 + coloc 后验比较 + 真实基因座 Mirror Regional Plots
# ==============================================================================
cis_pos[, label := sprintf("%s -> %s (nSNP=%d)", gene_symbol, gsub("BRCA_", "", outcome_label), nsnp)]
cis_pos[, sig_tag := fifelse(gene_symbol == "CXCR2", "LD-Confounded (PP.H3=1.0)",
                     fifelse(FDR < 0.05, "FDR < 0.05",
                     fifelse(gene_symbol %in% c("CXCL8", "PPARG"), "Plan A Core Axis (P<0.05)", "Nominal P < 0.05")))]
cis_pos[, domain_flag := fifelse(gene_symbol == "CXCR2", "LD-Confounded (Excluded)", "Primary")]
cis_pos[, anno_txt := sprintf("OR=%.3f [%.3f-%.3f], P=%.2g", OR, OR_LCI95, OR_UCI95, pval)]
cis_pos <- cis_pos[order(pval)]
cis_pos[, label := factor(label, levels = rev(label))]

p3a <- plot_forest_panel(
  cis_pos,
  "A. Strict cis-eQTLGen -> BCAC MR (8 Genes, 13 Pairs)",
  "Filtered to +-1 Mb cis-IVs (29 harmonized genes / 87 total tests); Steiger direction = TRUE for all",
  c(0.45, max(cis_pos$OR_UCI95) * 1.48),
  c("FDR < 0.05" = pal_npg[["red"]],
    "Plan A Core Axis (P<0.05)" = pal_npg[["green"]],
    "Nominal P < 0.05" = pal_npg[["blue"]],
    "LD-Confounded (PP.H3=1.0)" = "#94A3B8")
)

coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
coloc_dt[, pair := paste0(gene_symbol, " (", gsub("BRCA_", "", outcome_label), ")")]
coloc_long <- melt(coloc_dt, id.vars = "pair", measure.vars = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
                   variable.name = "Hypothesis", value.name = "PP")
coloc_dt <- coloc_dt[order(PP.H4)]
coloc_long[, pair := factor(pair, levels = coloc_dt$pair)]

p3b <- ggplot(coloc_long, aes(x = PP, y = pair, fill = Hypothesis)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("PP.H0"="#E2E8F0","PP.H1"="#CBD5E1","PP.H2"="#64748B","PP.H3"="#F97316","PP.H4"="#059669")) +
  labs(title = "B. Bayesian Colocalization (eQTLGen vs BCAC, +-500 kb)",
       subtitle = "IL1B ER- PP.H4=0.507 (H4/[H3+H4]=0.895); CXCR2 eliminated by PP.H3=0.9996",
       x = "Posterior Probability", y = NULL) +
  theme_journal() + theme(legend.position = "bottom")

# 3C: Regional Mirror Manhattan Comparison for 4 Representative Loci
reg_snps <- readRDS(file.path(cache_dir, "planA_coloc_plot_snps.rds"))
lead_sum <- fread(file.path(res_dir, "planA_coloc_lead_snps_summary.csv"))

plot_one_locus <- function(pair_key, title_lbl) {
  sub <- reg_snps[pair == pair_key]
  ls  <- lead_sum[pair == pair_key]
  sub[, pos_mb := pos / 1e6]
  top_df <- data.table(pos_mb = sub$pos_mb, val = sub$log10p_e, track = "eQTLGen -log10(P)", rsid = sub$rsid)
  bot_df <- data.table(pos_mb = sub$pos_mb, val = sub$log10p_o, track = "BCAC GWAS -log10(P)", rsid = sub$rsid)
  comb <- rbind(top_df, bot_df)
  comb[, track := factor(track, levels = c("eQTLGen -log10(P)", "BCAC GWAS -log10(P)"))]
  comb[, is_lead := (track == "eQTLGen -log10(P)" & rsid == ls$lead_eqtl_snp) | (track == "BCAC GWAS -log10(P)" & rsid == ls$lead_brca_snp)]

  ggplot(comb, aes(x = pos_mb, y = val)) +
    geom_point(data = comb[is_lead == FALSE], aes(color = track), size = 0.7, alpha = 0.45) +
    geom_point(data = comb[is_lead == TRUE], color = "#DC2626", shape = 18, size = 3.0) +
    ggrepel::geom_text_repel(data = comb[is_lead == TRUE], aes(label = rsid), size = 2.6, fontface = "bold", color = "#991B1B",
                             max.overlaps = 10, box.padding = 0.25) +
    facet_grid(track ~ ., scales = "free_y") +
    scale_color_manual(values = c("eQTLGen -log10(P)" = "#3C5488", "BCAC GWAS -log10(P)" = "#00A087"), guide = "none") +
    labs(title = title_lbl,
         subtitle = sprintf("Lead dist = %s kb (n=%d SNPs)", format(round(ls$bp_distance / 1000, 1), big.mark=","), ls$n_shared_snps),
         x = sprintf("Chr%s Position (Mb)", sub$chr[1]), y = "-log10(P)") +
    theme_journal(9.2)
}

p3c1 <- plot_one_locus("PPARG_ER_pos",  "C1. PPARG (ER+): Aligned Peaks")
p3c2 <- plot_one_locus("IL1B_ER_neg",   "C2. IL1B (ER-): PP.H4=0.507")
p3c3 <- plot_one_locus("CXCL8_ER_neg",  "C3. CXCL8 (ER-): Sub-thresh GWAS")
p3c4 <- plot_one_locus("CXCR2_Overall", "C4. CXCR2 (Excl.): PP.H3=1.0")

fig3_bot <- wrap_plots(p3c1, p3c2, p3c3, p3c4, ncol = 4)
fig3 <- (p3a | p3b) / fig3_bot + plot_layout(heights = c(1.0, 1.15))
save_pub(fig3, "Fig3_Pub_CisMR_Coloc_and_RegionalPeaks", w = 17.5, h = 12.0)

# ==============================================================================
# Fig 4: 方案 A 核心临床验证大图 (TCGA / METABRIC / 3x GEO + Meta-HR)
# ==============================================================================
tc  <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb  <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20 <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42 <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]

# A: PPARG Tumor vs Normal in TCGA & GSE42568
df_tn <- rbind(
  data.table(Cohort = "TCGA-BRCA\n(P=7.74e-57, FC=-3.23)", Group = "Normal (n=114)", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA\n(P=7.74e-57, FC=-3.23)", Group = "Tumor (n=1082)", Expr = as.numeric(tc$tumor_log2["PPARG", ])),
  data.table(Cohort = "GSE42568\n(P=1.29e-08, FC=-4.74)", Group = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal (n=17)", "Tumor (n=104)"),
             Expr = as.numeric(g42$expr_target["PPARG", ]))
)
df_tn[, G_simple := ifelse(grepl("Normal", Group), "Normal", "Tumor")]
p4a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple)) +
  geom_violin(alpha = 0.35, trim = FALSE) +
  geom_boxplot(width = 0.28, outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#38BDF8", "Tumor" = "#E64B35")) +
  labs(title = "A. PPARG Downregulation in Breast Cancer vs Normal",
       x = NULL, y = "Normalized log2 Expression") +
  theme_journal() + theme(legend.position = "none")

# B: PAM50 5-Subtype Distribution for PPARG & CXCL8 in TCGA
pam_tc <- rbind(
  data.table(Gene = "PPARG (KW P=5.09e-16)", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ])),
  data.table(Gene = "CXCL8 (KW P=9.85e-33)", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ]))
)[Subtype %in% c("Normal","LumA","LumB","Her2","Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal","LumA","LumB","Her2","Basal"))]

p4b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_boxplot(width = 0.55, outlier.size = 0.5, alpha = 0.88) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal"="#94A3B8","LumA"="#3B82F6","LumB"="#6366F1","Her2"="#F59E0B","Basal"="#E11D48")) +
  labs(title = "B. PAM50 Subtype Expression (TCGA, n=984)",
       subtitle = "PPARG highest in Normal/LumA; CXCL8 Basal > LumA (P<1e-16)",
       x = "PAM50 Subtype", y = "log2(RSEM + 1)") +
  theme_journal() + theme(legend.position = "none")

# C: TCGA ER+ PPARG OS KM
idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp], expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4c <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_tc_p), data = df_tc_p, pval = TRUE,
  legend.title = "PPARG", legend.labs = c("Low", "High"),
  palette = c("#E64B35", "#3C5488"), title = sprintf("C. TCGA ER+ OS: PPARG (n=%d)", nrow(df_tc_p)),
  xlab = "Months", ylab = "Overall Survival", ggtheme = theme_journal()
)$plot + labs(subtitle = "Cox HR = 0.820 (P=0.0135) | Log-rank P = 0.00084")

# D: METABRIC ER+ PPARG OS KM
idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp], expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4d <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_mb_p), data = df_mb_p, pval = TRUE,
  legend.title = "PPARG", legend.labs = c("Low", "High"),
  palette = c("#E64B35", "#3C5488"), title = sprintf("D. METABRIC ER+ OS: PPARG (n=%d)", nrow(df_mb_p)),
  xlab = "Months", ylab = "Overall Survival", ggtheme = theme_journal()
)$plot + labs(subtitle = "Continuous Cox HR = 0.894 (P=0.00638) | Median KM P = 0.081")

# E: GSE20685 CXCL8 OS KM
pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4e <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_g20_c), data = df_g20_c, pval = TRUE,
  legend.title = "CXCL8", legend.labs = c("Low", "High"),
  palette = c("#3C5488", "#E64B35"), title = sprintf("E. GSE20685 OS: CXCL8 (n=%d)", nrow(df_g20_c)),
  xlab = "Years", ylab = "Overall Survival", ggtheme = theme_journal()
)$plot + labs(subtitle = "Cox HR = 1.251 (P=0.00057) | MFS P = 0.00163")

# F: METABRIC ER- CXCL8 OS KM
idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern], expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4f <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_mb_c), data = df_mb_c, pval = TRUE,
  legend.title = "CXCL8", legend.labs = c("Low", "High"),
  palette = c("#3C5488", "#E64B35"), title = sprintf("F. METABRIC ER- OS: CXCL8 (n=%d)", nrow(df_mb_c)),
  xlab = "Months", ylab = "Overall Survival", ggtheme = theme_journal()
)$plot + labs(subtitle = "Log-rank P = 0.0158 | Multivar OS (Age+Stage) P = 0.000727")

# G: 5-Cohort Multi-Cohort Forest with Pooled Meta-HR Summary
surv_all <- fread(file.path(res_dir, "module4_survival_analysis.csv"))
mv_cox   <- fread(file.path(res_dir, "module4_multivariable_cox.csv"))
meta_hr  <- fread(file.path(res_dir, "planA_multicohort_meta_survival.csv"))

hr_sub <- rbind(
  surv_all[gene_symbol == "PPARG" & ((cohort == "TCGA_BRCA" & stratum %in% c("Overall","ER_pos") & endpoint == "OS") |
                                     (cohort == "METABRIC" & stratum %in% c("Overall","ER_pos") & endpoint == "OS") |
                                     (cohort %in% c("GSE42568","GSE20685","GSE7390") & endpoint == "OS")),
           .(gene_symbol, label = paste(cohort, stratum, endpoint), HR, HR_LCI95, HR_UCI95, pval = cox_pval, is_meta = FALSE)],
  meta_hr[analysis %in% c("PPARG_5Cohort_Overall_OS", "PPARG_2Cohort_ERpos_OS"),
          .(gene_symbol = "PPARG", label = paste0("POOLED ", analysis, " (I2=", I2_pct, "%)"),
            HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, is_meta = TRUE)],
  surv_all[gene_symbol == "CXCL8" & ((cohort == "METABRIC" & stratum == "Overall" & endpoint %in% c("OS","RFS")) |
                                     (cohort %in% c("GSE42568","GSE20685","GSE7390") & endpoint %in% c("OS","MFS","RFS")) |
                                     (cohort == "TCGA_BRCA" & stratum == "Overall" & endpoint == "PFS")),
           .(gene_symbol, label = paste(cohort, stratum, endpoint), HR, HR_LCI95, HR_UCI95, pval = cox_pval, is_meta = FALSE)],
  mv_cox[gene_symbol == "CXCL8" & cohort == "METABRIC" & stratum == "Overall",
         .(gene_symbol, label = "METABRIC Overall OS_Multivar(Age+Stage)", HR, HR_LCI95, HR_UCI95, pval, is_meta = FALSE)],
  meta_hr[analysis %in% c("CXCL8_5Cohort_Progression_RFS_MFS_PFS", "CXCL8_4Ext_RFS_MFS"),
          .(gene_symbol = "CXCL8", label = paste0("POOLED ", analysis, " (I2=", I2_pct, "%)"),
            HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, is_meta = TRUE)]
)
hr_sub[, tag := fifelse(is_meta, "Pooled Meta-HR", fifelse(pval < 0.05, "Cohort P < 0.05", "Cohort P >= 0.05"))]
hr_sub[, ord_lbl := factor(label, levels = rev(unique(label)))]

p4g <- ggplot(hr_sub, aes(x = HR, y = ord_lbl, color = tag, shape = is_meta)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.25, linewidth = 0.7) +
  geom_point(aes(size = is_meta)) +
  geom_text(aes(x = max(HR_UCI95, na.rm = TRUE) * 1.03, label = sprintf("HR=%.3f, P=%.2g", HR, pval)),
            hjust = 0, size = 2.8, color = "#1E293B", show.legend = FALSE) +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Pooled Meta-HR" = "#7C3AED", "Cohort P < 0.05" = "#E64B35", "Cohort P >= 0.05" = "#94A3B8")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_manual(values = c("FALSE" = 2.6, "TRUE" = 4.0), guide = "none") +
  coord_cartesian(xlim = c(0.42, max(hr_sub$HR_UCI95) * 1.38)) +
  labs(title = "G. Cross-Cohort Prognostic Meta-Analysis (5 Cohorts, n = 3,677)",
       subtitle = "PPARG ER+ Meta-HR = 0.878 (P=3.8e-4, I2=0%); CXCL8 Progression Meta-HR = 1.109 (P=0.0048)",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal() + theme(legend.position = "bottom")

rm(tc, mb, g20, g42); gc()
fig4_top <- (p4a | p4b)
fig4_mid <- (p4c | p4d) / (p4e | p4f)
fig4_all <- fig4_top / fig4_mid / p4g + plot_layout(heights = c(0.9, 1.5, 1.3))
save_pub(fig4_all, "Fig4_Pub_Clinical_DualAxis_5Cohorts", w = 15.8, h = 18.5)

# ==============================================================================
# Fig 5: HPA 蛋白 IHC 实拍 + 单细胞 nCPM + TCGA TME & 检查点 + Hallmark GSEA
# ==============================================================================
sc_hpa <- fread(file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
sc_keep <- c("Neutrophils", "Macrophages", "Fibroblasts", "Adipocytes",
             "Lymphatic endothelial cells", "Vascular endothelial cells",
             "Breast secretory cells", "Breast hormone-responsive cells",
             "T-cells", "B-cells", "Nk-cells")
sc_sub <- sc_hpa[gene_symbol %in% c("CXCL8", "PPARG") & cell_type %in% sc_keep]
sc_sub[, log10_ncpm := log10(nCPM + 1)]
sc_ord <- sc_sub[, .(mx = max(log10_ncpm)), by = cell_type][order(mx)]$cell_type
sc_sub[, cell_type := factor(cell_type, levels = sc_ord)]

p5a <- ggplot(sc_sub, aes(x = log10_ncpm, y = cell_type, fill = gene_symbol)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  scale_fill_manual(values = c("CXCL8" = pal_npg[["red"]], "PPARG" = pal_npg[["blue"]])) +
  labs(title = "A. HPA Single-Cell Type Expression (nCPM)",
       subtitle = "CXCL8 enriched in Neutrophils (5782) & Macrophages (2093); PPARG in Adipocytes (2599) & Endothelial (234)",
       x = "log10(nCPM + 1)", y = NULL, fill = "Gene") +
  theme_journal() + theme(legend.position = "bottom")

# 5B: Embed 3 official HPA IHC images cleanly via cowplot/ggplot raster
img_panel <- function(jpg_path, ttl, sub_t) {
  img <- jpeg::readJPEG(jpg_path)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_fixed(ratio = 1) +
    labs(title = ttl, subtitle = sub_t) +
    theme_void() +
    theme(aspect.ratio = 1,
          plot.title = element_text(face = "bold", size = 8.4, hjust = 0.5),
          plot.subtitle = element_text(size = 7.4, color = "#475569", hjust = 0.5),
          panel.border = element_rect(color = "#CBD5E1", fill = NA, linewidth = 0.6))
}
p5b1 <- img_panel(file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"), "B1. Normal PPARG", "Adipocyte Low | Duct ND")
p5b2 <- img_panel(file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"), "B2. BRCA PPARG IHC", "4 Med / 6 Low / 2 ND")
p5b3 <- img_panel(file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"), "B3. BRCA CXCL8 IHC", "12/12 ND (Secreted)")
p5b  <- (p5b1 | p5b2 | p5b3)

# 5C: TME & Immune Checkpoint Correlations
cor5 <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
tme_top <- cor5[((gene_symbol == "CXCL8" & stratum == "ER_neg") | (gene_symbol == "PPARG" & stratum == "ER_pos")) &
                  feature_type %in% c("TME_Cell_Score", "Checkpoint", "Genomic_Feature") & FDR < 0.05]
# Keep cleanest non-redundant set up to 12 per gene
tme_top <- tme_top[order(gene_symbol, -abs(rho))][, .SD[1:min(12, .N)], by = gene_symbol]
tme_top[, feat_clean := sub("^(quanTIseq|MCPcounter|ESTIMATE)__(.*)$", "\\2 [\\1]", feature)]
tme_top[, facet_lbl := ifelse(gene_symbol == "CXCL8", "C1. CXCL8 in ER- BRCA (n=237)", "C2. PPARG in ER+ BRCA (n=783)")]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feat_clean)]
tme_top <- tme_top[order(facet_lbl, rho)]
tme_top[, feat_ord := factor(feat_ord, levels = unique(feat_ord))]

p5c <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_segment(aes(x = 0, xend = rho, y = feat_ord, yend = feat_ord), linewidth = 0.8) +
  geom_point(size = 2.8) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("Positive" = pal_npg[["red"]], "Negative" = pal_npg[["blue"]])) +
  labs(title = "C. Immune Microenvironment, Checkpoints & Genomic Instability (FDR < 0.05)",
       x = "Spearman rho", y = NULL, color = NULL) +
  theme_journal() + theme(legend.position = "bottom")

# 5D: Hallmark GSEA
gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene == "CXCL8" & stratum == "ER_neg") | (target_gene == "PPARG" & stratum == "ER_pos")) & p.adjust < 0.05]
gs_top <- gs_top[order(target_gene, -abs(NES))][, .SD[1:min(10, .N)], by = target_gene]
gs_top[, clean_id := gsub("^HALLMARK_", "", ID)]
gs_top[, facet_lbl := ifelse(target_gene == "CXCL8", "D1. CXCL8 High vs Low (ER- BRCA)", "D2. PPARG High vs Low (ER+ BRCA)")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]
gs_top[, id_ord := factor(id_ord, levels = id_ord)]

p5d <- ggplot(gs_top, aes(x = NES, y = id_ord, fill = -log10(p.adjust))) +
  geom_col(width = 0.72) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_fill_gradient(low = "#93C5FD", high = "#1E3A8A") +
  labs(title = "D. Single-Gene Hallmark GSEA (FDR < 0.05)",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = "-log10(FDR)") +
  theme_journal() + theme(legend.position = "bottom")

fig5 <- (p5a | p5b) / p5c / p5d + plot_layout(heights = c(0.95, 1.0, 1.0))
save_pub(fig5, "Fig5_Pub_HPA_TME_SingleCell_GSEA", w = 16.5, h = 14.5)

# ==============================================================================
# Fig 6: GDSC2 药敏与亚型去混杂散点验证 (47 株乳腺癌细胞系)
# ==============================================================================
sub_d <- fread(file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
d_dedup <- sub_d[order(gene_symbol, all_spearman_pval)][, .SD[1], by = .(gene_symbol, DRUG_NAME)]
d_top <- d_dedup[all_spearman_pval < 0.05][order(gene_symbol, all_spearman_pval)][, .SD[1:min(10, .N)], by = gene_symbol]
d_top[, dir := fifelse(subtype_adjusted_pval < 0.05,
               fifelse(all_spearman_rho < 0, "Sensitive + Subtype-Adj P<0.05", "Resistant + Subtype-Adj P<0.05"),
               fifelse(all_spearman_rho < 0, "Sensitive (Overall P<0.05)", "Resistant (Overall P<0.05)"))]
d_top[, drug_ord := paste0(gene_symbol, "___", DRUG_NAME)]
d_top <- d_top[order(gene_symbol, -all_spearman_pval)]
d_top[, drug_ord := factor(drug_ord, levels = drug_ord)]

p6a <- ggplot(d_top, aes(x = all_spearman_rho, y = drug_ord, color = dir, size = -log10(all_spearman_pval))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_point() +
  facet_wrap(~ gene_symbol, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c(
    "Sensitive + Subtype-Adj P<0.05" = "#059669",
    "Sensitive (Overall P<0.05)"     = "#6EE7B7",
    "Resistant + Subtype-Adj P<0.05" = "#DC2626",
    "Resistant (Overall P<0.05)"     = "#FCA5A5"
  )) +
  labs(title = "A. GDSC2 Drug Sensitivity (n=47 BRCA Cell Lines) & Subtype-Adjusted Robustness",
       subtitle = "Dark colors survive OLS adjustment for ER+/HER2+/TNBC subtype (Fulvestrant/Vorinostat/POMHEX adj P<0.05)",
       x = "Spearman rho (Expression vs LN_IC50)", y = NULL, color = "Response & Subtype Test", size = "-log10(P)") +
  theme_journal() + theme(legend.position = "bottom")

# 6B & 6C: Scatter Plots Colored by Cell Line Subtype
gb_mat <- readRDS(file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))
sc_das <- gb_mat[DRUG_NAME == "Dasatinib"]
p6b <- ggplot(sc_das, aes(x = CXCL8, y = LN_IC50, color = Subtype_Group)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(data = sc_das[Subtype_Group == "TNBC"], method = "lm", se = FALSE, color = "#E11D48", linewidth = 0.85) +
  scale_color_manual(values = c("ER_pos" = "#3B82F6", "HER2_pos" = "#F59E0B", "TNBC" = "#E11D48")) +
  labs(title = "B. CXCL8 vs Dasatinib Sensitivity",
       subtitle = "All (n=46): rho=-0.466, P=0.0011 | Within TNBC (n=21): rho=-0.417, P=0.060",
       x = "CXCL8 Expression (log2 TPM+1)", y = "Dasatinib LN_IC50", color = "Subtype") +
  theme_journal() + theme(legend.position = "bottom")

sc_fulv <- gb_mat[DRUG_ID == 1200] # Fulvestrant n=46 (adj_P=0.0096)
p6c <- ggplot(sc_fulv, aes(x = PPARG, y = LN_IC50, color = Subtype_Group)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(data = sc_fulv[Subtype_Group == "ER_pos"], method = "lm", se = FALSE, color = "#3B82F6", linewidth = 0.85) +
  scale_color_manual(values = c("ER_pos" = "#3B82F6", "HER2_pos" = "#F59E0B", "TNBC" = "#E11D48")) +
  labs(title = "C. PPARG vs Fulvestrant Resistance",
       subtitle = "All (n=46): rho=+0.375, P=0.011 | Subtype-Adj P=0.0096 | ER+ (n=14): rho=+0.503, P=0.069",
       x = "PPARG Expression (log2 TPM+1)", y = "Fulvestrant LN_IC50", color = "Subtype") +
  theme_journal() + theme(legend.position = "bottom")

save_pub(p6a / (p6b | p6c) + plot_layout(heights = c(1.1, 1.0)), "Fig6_Pub_Pharmacogenomics_Subtyped", w = 14.5, h = 10.5)

# ==============================================================================
# Fig S1: 上游核心因果菌群 & 代谢物 MR 多方法散点图 + Leave-One-Out 敏感性组合图
# ==============================================================================
single_dt <- fread(file.path(res_dir, "planA_mr_single_snp_wald.csv"))
sens_dt   <- fread(file.path(res_dir, "planA_mr_sensitivity_summary.csv"))
loo_dt    <- fread(file.path(res_dir, "planA_mr_leaveoneout.csv"))

target_s1 <- list(
  list(eid = "ebi-a-GCST90016960", oid = "ieu-a-1127", lbl = "Adlercreutzia -> ER+"),
  list(eid = "ebi-a-GCST90016987", oid = "ieu-a-1128", lbl = "Desulfovibrio -> ER-"),
  list(eid = "ebi-a-GCST90017039", oid = "ieu-a-1126", lbl = "Parabacteroides -> Overall"),
  list(eid = "ebi-a-GCST90016951", oid = "ieu-a-1127", lbl = "Ruminococcaceae -> ER+"),
  list(eid = "ebi-a-GCST90016989", oid = "ieu-a-1128", lbl = "Dorea -> ER-"),
  list(eid = "met-d-Acetate",      oid = "ieu-a-1126", lbl = "Acetate -> Overall")
)

sc_plots  <- list()
loo_plots <- list()
for (tk in target_s1) {
  sd <- single_dt[id.exposure == tk$eid & id.outcome == tk$oid]
  # Normalize SNP orientation so beta.exposure > 0
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
    scale_color_manual(values = c("Inverse variance weighted" = "#E64B35", "Weighted median" = "#00A087",
                                  "MR Egger" = "#3C5488", "Weighted mode" = "#F59E0B")) +
    labs(title = tk$lbl, x = "SNP Effect on Exposure", y = "SNP Effect on BRCA", color = "MR Method") +
    theme_journal(9.0) + theme(legend.position = "none")
  sc_plots[[tk$lbl]] <- p_sc

  ld <- loo_dt[id.exposure == tk$eid & id.outcome == tk$oid]
  ld[, is_all := (SNP == "All")]
  ld <- ld[order(is_all, -b)]
  ld[, snp_ord := factor(SNP, levels = unique(rev(SNP)))]
  p_lo <- ggplot(ld, aes(x = b, y = snp_ord, color = is_all)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
    geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), orientation = "y", width = 0.25, linewidth = 0.55) +
    geom_point(size = 1.9) +
    scale_color_manual(values = c("FALSE" = "#475569", "TRUE" = "#DC2626"), guide = "none") +
    labs(title = paste("LOO:", tk$lbl), x = "Leave-one-out IVW beta", y = NULL) +
    theme_journal(9.0)
  loo_plots[[tk$lbl]] <- p_lo
}

fig_s1_top <- wrap_plots(sc_plots, ncol = 3) +
  plot_annotation(title = "Supplementary Figure S1A: MR Multi-Method Scatter Plots for Top Causal Exposures")
# Extract single legend cleanly from first pair
fig_s1_bot <- wrap_plots(loo_plots, ncol = 3) +
  plot_annotation(title = "Supplementary Figure S1B: Leave-One-Out Sensitivity Analysis (All Leave-One-Out Directions Consistent)")
fig_s1 <- wrap_plots(c(sc_plots, loo_plots), ncol = 3) +
  plot_annotation(title = "Supplementary Figure S1: MR Multi-Method Concordance (Top 6 Panels) & Leave-One-Out Stability (Bottom 6 Panels)",
                  subtitle = "Lines: Red=IVW, Green=Weighted Median, Blue=MR-Egger, Amber=Weighted Mode; Bottom red points show overall IVW beta")
save_pub(fig_s1, "FigS1_Pub_Upstream_MR_Sensitivity", w = 15.5, h = 14.5)

cat("\nALL Plan A Publication Figures Exported Successfully to:", pub_dir, "\n")
