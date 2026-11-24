# 21_generate_final_planA_figures.R
# 方案 A 抗审稿完整投稿级高清组图总控脚本 (Fig 1 - Fig 6 + Fig S1 - Fig S2)
# 完整整合: 全局质控漏斗 + DMP/FR02跨队列 + IL-8 pQTL蛋白层 + 脂肪校正/Nomogram/异质性分解 + 真实肿瘤单细胞 + 临床治疗队列 + GDSC2

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(survminer)
  library(rms)
  library(jpeg)
  library(png)
  library(grid)
  library(scales)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
bulk_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
sc_dir    <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/scRNA_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
pub_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/Publication_PlanA")
ihc_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/HPA_IHC_Images")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

pal_npg <- c(
  red   = "#E64B35",
  blue  = "#3C5488",
  cyan  = "#4DBBD5",
  green = "#00A087",
  coral = "#F39B7F",
  slate = "#8491B4",
  amber = "#D97706"
)

theme_journal <- function(base_size = 10.0) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#F1F5F9", linewidth = 0.4),
      panel.border = element_rect(color = "#334155", linewidth = 0.65),
      strip.background = element_rect(fill = "#F8FAFC", color = "#CBD5E1"),
      strip.text = element_text(face = "bold", size = base_size - 0.5, color = "#0F172A"),
      plot.title = element_text(face = "bold", size = base_size + 1.2, color = "#0F172A"),
      plot.subtitle = element_text(size = base_size - 1.3, color = "#475569"),
      axis.title = element_text(face = "bold", size = base_size - 0.6, color = "#1E293B"),
      axis.text = element_text(color = "#334155", size = base_size - 1.3),
      legend.title = element_text(face = "bold", size = base_size - 1.0),
      legend.text = element_text(size = base_size - 1.3),
      plot.tag = element_text(face = "bold", size = base_size + 2.5, color = "#0F172A")
    )
}

save_pub <- function(p, name, w = 15, h = 8) {
  ggsave(file.path(pub_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(pub_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  cat("Exported:", name, "\n")
}

plot_forest_panel <- function(dt, title_txt, sub_txt, xlim_vals, color_map, text_x_mult = 1.04) {
  ggplot(dt, aes(x = OR, y = label, color = sig_tag)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.6) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.72) +
    geom_point(aes(shape = domain_flag), size = 2.8) +
    geom_text(aes(x = max(OR_UCI95, na.rm = TRUE) * text_x_mult, label = anno_txt),
              hjust = 0, size = 2.75, color = "#1E293B", show.legend = FALSE) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = c("Primary" = 16, "Warning (2 SNPs)" = 17, "LD-Confounded (Excluded)" = 4), guide = "none") +
    coord_cartesian(xlim = xlim_vals) +
    labs(title = title_txt, subtitle = sub_txt, x = "Odds Ratio (95% CI)", y = NULL, color = "Evidence Tier") +
    theme_journal() +
    theme(legend.position = "bottom")
}

# ==============================================================================
# FIG 1: 全局质控漏斗 (1A) + 上游菌群/跨队列Meta (1B) + 血浆代谢物 MR (1C)
# ==============================================================================
boxes <- data.table(
  id = c("s1", "qc1", "s2", "qc2", "s3", "qc3", "qc4", "s4_erneg", "s4_erpos"),
  xmin = c( 2,  2, 35, 35, 68, 68, 68,  8, 54),
  xmax = c(30, 30, 63, 63, 98, 98, 98, 46, 94),
  ymin = c(44, 30, 44, 30, 44, 31, 18,  2,  2),
  ymax = c(56, 40, 56, 40, 56, 40, 27, 13, 13),
  fill = c("#EFF6FF","#FEF2F2","#F0FDF4","#FFFBEB","#F5F3FF","#FEF2F2","#FFF7ED","#FFF1F2","#EFF6FF"),
  col  = c("#1D4ED8","#DC2626","#15803D","#D97706","#6D28D9","#DC2626","#EA580C","#BE123C","#1E40AF"),
  lbl  = c(
    "1. Microbiome Two-Sample MR\n211 MiBioGen (196 valid) -> 64 aligned (192 tests)\n13 nominal pairs (7 taxa) -> 6 validated taxa",
    "QC Gate 1: Reverse MR\nExcluded Marvinbryantia (rev IVW P=0.0216)\nDMP_2022(84)+FR02(36) Meta: Dorea/Desulfovibrio P<0.05",
    "2. Plasma Metabolome MR\ngutMGene bioactive intermediates -> 47 OpenGWAS\n141 tests -> 5 priority metabolites (Acetate FDR<0.05)",
    "QC Gate 2: Reverse MR & FinnGen\nAll 15 reverse-MR IVW P>0.05; Linoleic acid\nER- P=0.031 (Step-1 taxa->met P<0.05)",
    "3. Host Causal Targets (eQTL + pQTL)\n61 gutMGene genes -> 41 eQTLGen (123 tests)\n8 cis-positive genes (13 pairs, Steiger TRUE)",
    "QC Gate 3: Strict +-1 Mb cis-IV Filter\nExcluded 8 trans-only genes (FOXP3 P=2.67e-8)\nRetained 29 cis-harmonized genes (87 tests)",
    "QC Gate 4: Bayesian Coloc + Plasma pQTL\nExcluded CXCR2 (PP.H3=0.9996; IL1B H4=0.507)\nSCALLOP+INTERVAL cis-pQTL -> IL-8 ER- (P=0.036/0.033)",
    "Plan A ER- Risk Axis: CXCL8 (+IL1B/FFAR2)\nDorea/Parabacteroides -| CXCL8 <- Acetate | 3-GEO Meta I2=0% P=2.1e-6\nTISCH2 Mono/Macro (31%); GSE25066 DRFS P=0.003; I-SPY2 pCR Pint=0.0056",
    "Plan A ER+ Protective Axis: PPARG\nAdlercreutzia -> Equol & Linoleic->10-keto-12Z-octadecenoic->PPARG (PMID:25749343)\nTCGA OS indep. of ADIPOQ+PLIN1+Stromal (P=0.032) & Nomogram (P=0.028)\nTISCH2 Adipocyte-Free Endothelial/Malignant; METABRIC HT=YES OS P=0.0065"
  )
)

arrows_df <- data.table(
  x1 = c(16, 49, 83, 83, 16, 22, 49, 49, 83, 83),
  y1 = c(44, 44, 44, 31, 30, 30, 30, 30, 18, 18),
  x2 = c(16, 49, 83, 83, 27, 68, 27, 74, 27, 74),
  y2 = c(40, 40, 40, 27, 13, 13, 13, 13, 13, 13)
)

p1a <- ggplot() +
  geom_segment(data = arrows_df, aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "#475569", linewidth = 0.75, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = col),
            linewidth = 0.75, show.legend = FALSE) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = lbl, color = col),
            size = 2.65, fontface = "bold", lineheight = 1.12, show.legend = FALSE) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 58), expand = FALSE) +
  labs(title = "A. Multi-Stage Causal Discovery & Quality-Control Funnel (Reverse-MR, cis-IV, Colocalization & Subtype Decoupling)",
       subtitle = "Sequential elimination of reverse causation (Marvinbryantia), trans-pleiotropy (FOXP3), LD confounding (CXCR2), and bulk adiposity artifacts") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 11.2, color = "#0F172A", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 8.8, color = "#475569", margin = margin(b = 6))
  )

mib_sig <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))
mib_sig <- mib_sig[gwas_name != "Marvinbryantia"]
mib_sig[, label := sprintf("%s %s -> %s (nSNP=%d)", gwas_rank, gwas_name, gsub("BRCA_", "", outcome_label), nsnp)]
mib_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", fifelse(FDR < 0.10, "FDR < 0.10", "Nominal P < 0.05"))]
mib_sig[, domain_flag := fifelse(nsnp <= 2, "Warning (2 SNPs)", "Primary")]
mib_sig[, anno_txt := sprintf("OR=%.2f [%.2f-%.2f], P=%.2g", OR, OR_LCI95, OR_UCI95, pval)]
mib_sig <- mib_sig[order(pval)]
mib_sig[, label := factor(label, levels = rev(label))]

p1b <- plot_forest_panel(
  mib_sig,
  "B. Causal Gut Microbiota -> BCAC (MiBioGen Validated)",
  "6 taxa after Marvinbryantia removal; Dorea (ER- Meta P=0.021) & Desulfovibrio (P=0.022) replicate",
  c(0.62, max(mib_sig$OR_UCI95) * 1.48),
  c("FDR < 0.05" = pal_npg[["red"]], "FDR < 0.10" = pal_npg[["amber"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

met_sig <- fread(file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv"))
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

p1c <- plot_forest_panel(
  met_sig,
  "C. Causal Plasma Metabolites -> BCAC",
  "5 prioritized gutMGene metabolites (all 15 reverse-MR IVW P>0.05; Acetate FDR<0.05)",
  c(0.45, max(met_sig$OR_UCI95) * 1.50),
  c("FDR < 0.05" = pal_npg[["red"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

fig1_all <- p1a / (p1b | p1c) + plot_layout(heights = c(1.05, 1.0))
save_pub(fig1_all, "Fig1_Pub_Design_QC_and_Upstream_MR", w = 17.5, h = 11.8)
rm(fig1_all, p1a, p1b, p1c); gc()

# ==============================================================================
# FIG 2: 方案 A 双主轴多组学因果网络 (含 DMP/FR02 与 IL-8 cis-pQTL 标签升级)
# ==============================================================================
edges_dt <- fread(file.path(res_dir, "module3_candidate_gutmgene_edges.csv"))
cis_pos  <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_g    <- unique(cis_pos$gene_symbol)

sub_e <- edges_dt[gene_symbol %in% pos_g & !is.na(metabolite)]
map_microbe <- function(m) {
  if (grepl("Parabacteroides", m, ignore.case = TRUE)) return("Parabacteroides (Prot.)")
  if (grepl("Dorea", m, ignore.case = TRUE)) return("Dorea (ER- Meta P=0.021)")
  if (grepl("Ruminococc", m, ignore.case = TRUE)) return("Ruminococcaceae (Risk)")
  if (grepl("Blautia", m, ignore.case = TRUE)) return("Blautia (ER- Risk, 2 SNPs)")
  if (grepl("MR_Positive_Metabolite", m)) return("Blood Acetate (Risk, FDR<0.05)")
  return(NA_character_)
}
sub_e[, m_grp := sapply(microbe, map_microbe)]
mm_links <- unique(sub_e[!is.na(m_grp), .(from = m_grp, to = metabolite)])
mm_links <- unique(rbind(
  mm_links,
  data.table(from = "Adlercreutzia (ER+ Meta P=0.0056)", to = "Equol"),
  data.table(from = "Blood Acetate (Risk, FDR<0.05)", to = "Acetate"),
  data.table(from = "Blood Linoleic Acid (ER- MR, P=0.031)", to = "10-keto-12Z-octadecenoic acid")
))
mg_links <- unique(sub_e[, .(from = metabolite, to = gene_symbol, alteration)])
mg_links <- mg_links[, .(alteration = ifelse(length(unique(alteration)) > 1, "mixed", alteration[1])), by = .(from, to)]
mg_links <- unique(rbind(
  mg_links,
  data.table(from = "10-keto-12Z-octadecenoic acid", to = "PPARG", alteration = "activation")
))

left_nodes  <- c("Parabacteroides (Prot.)", "Dorea (ER- Meta P=0.021)", "Blautia (ER- Risk, 2 SNPs)",
                 "Blood Acetate (Risk, FDR<0.05)", "Ruminococcaceae (Risk)",
                 "Blood Linoleic Acid (ER- MR, P=0.031)", "Adlercreutzia (ER+ Meta P=0.0056)")
mid_nodes   <- c("Butyrate", "Acetate", "Succinate", "3-Indolepropionic acid",
                 "Trimethylamine oxide", "10-keto-12Z-octadecenoic acid", "Equol")
right_nodes <- c("CXCL8 [ER- eQTL + 2x pQTL P<0.04]", "IL1B [Coloc PP.H4=0.507]", "FFAR2 [GTEx Breast]",
                 "TNF [Secondary]", "HDAC11 [Secondary]", "IL12A [Secondary]",
                 "PPARG [ER+ eQTL + Adipose-Adj OS]", "CXCR2 [Excluded: PP.H3=1.0]")

rename_g <- c(
  "CXCL8"  = "CXCL8 [ER- eQTL + 2x pQTL P<0.04]",
  "IL1B"   = "IL1B [Coloc PP.H4=0.507]",
  "FFAR2"  = "FFAR2 [GTEx Breast]",
  "TNF"    = "TNF [Secondary]",
  "HDAC11" = "HDAC11 [Secondary]",
  "IL12A"  = "IL12A [Secondary]",
  "PPARG"  = "PPARG [ER+ eQTL + Adipose-Adj OS]",
  "CXCR2"  = "CXCR2 [Excluded: PP.H3=1.0]"
)
mg_links[, to := rename_g[to]]
mm_links <- mm_links[from %in% left_nodes & to %in% mid_nodes]
mg_links <- mg_links[from %in% mid_nodes & to %in% right_nodes]

node_pos <- rbindlist(list(
  data.table(node = left_nodes,  col = 1, y = seq(9.5, 1.5, length.out = length(left_nodes))),
  data.table(node = mid_nodes,   col = 2, y = seq(9.5, 1.5, length.out = length(mid_nodes))),
  data.table(node = right_nodes, col = 3, y = seq(9.5, 1.5, length.out = length(right_nodes)))
))
e1 <- merge(mm_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e1 <- merge(e1, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to"); e1[, edge_type := "produces"]
e2 <- merge(mg_links, node_pos[, .(from = node, x1 = col, y1 = y)], by = "from")
e2 <- merge(e2, node_pos[, .(to = node, x2 = col, y2 = y)], by = "to"); e2[, edge_type := alteration]
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
             size = 3.1, fontface = "bold", label.padding = unit(0.32, "lines"), linewidth = 0.6) +
  scale_fill_identity() +
  scale_color_manual(
    values = c("produces" = "#64748B", "activation" = "#059669", "inhibition" = "#DC2626", "mixed" = "#7C3AED",
               "#BE123C" = "#BE123C", "#1D4ED8" = "#1D4ED8", "#94A3B8" = "#94A3B8", "#0F172A" = "#0F172A"),
    breaks = c("produces", "activation", "inhibition", "mixed"), name = "gutMGene Edge"
  ) +
  scale_linetype_manual(
    values = c("produces" = "solid", "activation" = "solid", "inhibition" = "dashed", "mixed" = "dotted"),
    breaks = c("produces", "activation", "inhibition", "mixed"), name = "gutMGene Edge"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("Causal Microbe / Metabolite Arm (MR + Cross-Cohort)", "gutMGene Bioactive Mediators", "Strict cis-eQTL & Plasma pQTL Targets"),
                     limits = c(0.45, 3.58)) +
  ylim(0.8, 10.2) +
  labs(title = "Plan A Multi-Omics Causal Architecture: Subtype-Divergent Microbiota-Metabolite-Host Axes",
       subtitle = "ER- Axis (Rose): Dorea/Parabacteroides -> Butyrate/Acetate -| CXCL8 (mRNA + IL-8 pQTL); ER+ Axis (Blue): Linoleic acid -> 10-keto-12Z-octadecenoic -> PPARG",
       x = NULL, y = NULL) +
  theme_journal() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
save_pub(p2, "Fig2_Pub_PlanA_DualAxis_Network", w = 14.8, h = 7.6)
rm(p2, all_seg, e1, e2); gc()

# ==============================================================================
# FIG 3: cis-eQTL MR (3A) + IL-8 血浆 pQTL MR (3B) + Coloc (3C) + 基因座镜像图 (3D)
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
  "Filtered to +-1 Mb cis-IVs (29 harmonized genes / 87 total tests); Steiger direction = TRUE",
  c(0.45, max(cis_pos$OR_UCI95) * 1.48),
  c("FDR < 0.05" = pal_npg[["red"]], "Plan A Core Axis (P<0.05)" = pal_npg[["green"]],
    "Nominal P < 0.05" = pal_npg[["blue"]], "LD-Confounded (PP.H3=1.0)" = "#94A3B8")
)

# 3B: IL-8 Plasma pQTL -> BCAC MR across 5 cohorts (cis vs All IVs)
pqtl_dt <- fread(file.path(res_dir, "planA_pqtl_il8_mr.csv"))
pq_sub <- pqtl_dt[method %in% c("Inverse variance weighted", "Wald ratio") &
                  outcome_label %in% c("BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall")]
# Keep unique cohort + iv_scope (SCALLOP 5e-8 Cis == 1e-5 Cis, keep 5e-8 Cis to avoid duplicate)
pq_sub <- pq_sub[!(id.exposure == "ebi-a-GCST90011994" & thresh == "1e-5" & iv_scope == "Cis_1Mb_CXCL8")]
pq_sub[, cohort_lbl := fifelse(id.exposure == "ebi-a-GCST90011994" & iv_scope == "Cis_1Mb_CXCL8", "SCALLOP 2020 (n=21,758) cis-pQTL [rs7655660]",
                       fifelse(id.exposure == "ebi-a-GCST90011994", "SCALLOP 2020 (n=21,758) All 31 IVs",
                       fifelse(id.exposure == "prot-a-749" & iv_scope == "Cis_1Mb_CXCL8", "INTERVAL 2018 (n=3,301) cis-pQTL [rs74327761]",
                       fifelse(id.exposure == "prot-a-749", "INTERVAL 2018 (n=3,301) All 24 IVs",
                       fifelse(id.exposure == "ebi-a-GCST004445", "YFS/FINRISK 2016 (n=3,526) All 13 IVs",
                       fifelse(id.exposure == "prot-c-3447_64_2", "KORA 2019 (Suhre) All 7 IVs",
                               "LBC 2020 (n=982) All 10 IVs"))))))]
pq_sub[, tier := fifelse(pval < 0.05, "Significant (P < 0.05)", "Null (P >= 0.05)")]
pq_sub[, out_clean := factor(outcome_label, levels = c("BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall"),
                             labels = c("ER- Breast Cancer", "ER+ Breast Cancer", "Overall Breast Cancer"))]
ord_pq <- c("SCALLOP 2020 (n=21,758) cis-pQTL [rs7655660]", "INTERVAL 2018 (n=3,301) cis-pQTL [rs74327761]",
            "SCALLOP 2020 (n=21,758) All 31 IVs", "INTERVAL 2018 (n=3,301) All 24 IVs",
            "YFS/FINRISK 2016 (n=3,526) All 13 IVs", "LBC 2020 (n=982) All 10 IVs", "KORA 2019 (Suhre) All 7 IVs")
pq_sub[, cohort_lbl := factor(cohort_lbl, levels = rev(ord_pq))]

p3b <- ggplot(pq_sub, aes(x = OR, y = cohort_lbl, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.55) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.7) +
  geom_point(size = 2.7) +
  geom_text(aes(x = pmin(1.76, OR_UCI95 + 0.02), label = ifelse(pval < 0.05, sprintf("%.3f (P=%.3f)", OR, pval), "")),
            hjust = 0, size = 2.55, fontface = "bold", color = "#BE123C", show.legend = FALSE) +
  facet_wrap(~ out_clean, ncol = 3) +
  scale_color_manual(values = c("Significant (P < 0.05)" = pal_npg[["red"]], "Null (P >= 0.05)" = "#94A3B8")) +
  coord_cartesian(xlim = c(0.84, 1.98)) +
  labs(title = "B. Protein-Layer Validation: Plasma CXCL8/IL-8 pQTL -> BCAC MR (5 Cohorts)",
       subtitle = "Both SCALLOP (P=0.0363) & INTERVAL (P=0.0328) cis-pQTLs replicate ER--specific causality (PPARG: 0 OpenGWAS plasma pQTLs)",
       x = "Odds Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal() + theme(legend.position = "bottom")

# 3C: Colocalization including SCALLOP IL-8 pQTL row
coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
coloc_dt[, pair := paste0(gene_symbol, " (", gsub("BRCA_", "", outcome_label), ")")]
pq_coloc <- fread(file.path(res_dir, "planA_pqtl_il8_coloc.csv"))
coloc_comb <- rbind(
  coloc_dt[, .(pair, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)],
  data.table(pair = "IL-8 pQTL (SCALLOP ER-)", PP.H0 = pq_coloc$PP.H0.abf[1], PP.H1 = pq_coloc$PP.H1.abf[1],
             PP.H2 = pq_coloc$PP.H2.abf[1], PP.H3 = pq_coloc$PP.H3.abf[1], PP.H4 = pq_coloc$PP.H4.abf[1])
)
coloc_comb <- coloc_comb[order(PP.H4)]
coloc_long <- melt(coloc_comb, id.vars = "pair", measure.vars = c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),
                   variable.name = "Hypothesis", value.name = "PP")
coloc_long[, pair := factor(pair, levels = coloc_comb$pair)]

p3c <- ggplot(coloc_long, aes(x = PP, y = pair, fill = Hypothesis)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("PP.H0"="#E2E8F0","PP.H1"="#CBD5E1","PP.H2"="#64748B","PP.H3"="#F97316","PP.H4"="#059669")) +
  labs(title = "C. Bayesian Colocalization (+-500 kb)",
       subtitle = "IL1B ER- PP.H4=0.507 (H4/[H3+H4]=0.895); CXCR2 eliminated by PP.H3=0.9996",
       x = "Posterior Probability", y = NULL) +
  theme_journal() + theme(legend.position = "bottom")

# 3D: Regional Mirror Plots
reg_snps <- readRDS(file.path(cache_dir, "planA_coloc_plot_snps.rds"))
lead_sum <- fread(file.path(res_dir, "planA_coloc_lead_snps_summary.csv"))

plot_one_locus <- function(pair_key, title_lbl) {
  sub <- reg_snps[pair == pair_key]
  ls  <- lead_sum[pair == pair_key]
  sub[, pos_mb := pos / 1e6]
  comb <- rbind(
    data.table(pos_mb = sub$pos_mb, val = sub$log10p_e, track = "eQTLGen -log10(P)", rsid = sub$rsid),
    data.table(pos_mb = sub$pos_mb, val = sub$log10p_o, track = "BCAC GWAS -log10(P)", rsid = sub$rsid)
  )
  comb[, track := factor(track, levels = c("eQTLGen -log10(P)", "BCAC GWAS -log10(P)"))]
  comb[, is_lead := (track == "eQTLGen -log10(P)" & rsid == ls$lead_eqtl_snp) | (track == "BCAC GWAS -log10(P)" & rsid == ls$lead_brca_snp)]

  ggplot(comb, aes(x = pos_mb, y = val)) +
    geom_point(data = comb[is_lead == FALSE], aes(color = track), size = 0.65, alpha = 0.42) +
    geom_point(data = comb[is_lead == TRUE], color = "#DC2626", shape = 18, size = 2.9) +
    ggrepel::geom_text_repel(data = comb[is_lead == TRUE], aes(label = rsid), size = 2.5, fontface = "bold", color = "#991B1B",
                             max.overlaps = 10, box.padding = 0.45, point.padding = 0.35, min.segment.length = 0) +
    facet_grid(track ~ ., scales = "free_y") +
    scale_color_manual(values = c("eQTLGen -log10(P)" = "#3C5488", "BCAC GWAS -log10(P)" = "#00A087"), guide = "none") +
    labs(title = title_lbl,
         subtitle = sprintf("Lead dist = %s kb (n=%d SNPs)", format(round(ls$bp_distance / 1000, 1), big.mark=","), ls$n_shared_snps),
         x = sprintf("Chr%s Position (Mb)", sub$chr[1]), y = "-log10(P)") +
    theme_journal(9.0)
}

p3d <- wrap_plots(
  plot_one_locus("PPARG_ER_pos",  "D1. PPARG (ER+): Aligned Peaks"),
  plot_one_locus("IL1B_ER_neg",   "D2. IL1B (ER-): PP.H4=0.507"),
  plot_one_locus("CXCL8_ER_neg",  "D3. CXCL8 (ER-): Sub-thresh GWAS"),
  plot_one_locus("CXCR2_Overall", "D4. CXCR2 (Excl.): PP.H3=1.0"),
  ncol = 4
)
fig3 <- (p3a | p3c) / p3b / p3d + plot_layout(heights = c(1.05, 0.85, 1.05))
save_pub(fig3, "Fig3_Pub_CisMR_pQTL_Coloc_and_Peaks", w = 17.5, h = 15.0)
rm(fig3, p3a, p3b, p3c, p3d, reg_snps); gc()

# ==============================================================================
# FIG 4: 多队列生存 + 异质性分解 + PPARG 脂肪校正多因素Cox + C-index/timeROC + 校准曲线
# ==============================================================================
tc  <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb  <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20 <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42 <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]

df_tn <- rbind(
  data.table(Cohort = "TCGA-BRCA\n(P=7.74e-57, FC=-3.23)", G_simple = "Normal", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA\n(P=7.74e-57, FC=-3.23)", G_simple = "Tumor",  Expr = as.numeric(tc$tumor_log2["PPARG", ])),
  data.table(Cohort = "GSE42568\n(P=1.29e-08, FC=-4.74)", G_simple = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal", "Tumor"),
             Expr = as.numeric(g42$expr_target["PPARG", ]))
)
p4a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple)) +
  geom_violin(alpha = 0.35, trim = FALSE) +
  geom_boxplot(width = 0.28, outlier.size = 0.45, alpha = 0.85) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#38BDF8", "Tumor" = "#E64B35")) +
  labs(title = "A. PPARG Downregulation in Tumor vs Normal", x = NULL, y = "Normalized log2 Expr") +
  theme_journal(9.5) + theme(legend.position = "none")

pam_tc <- rbind(
  data.table(Gene = "PPARG (KW P=5.09e-16)", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ])),
  data.table(Gene = "CXCL8 (KW P=9.85e-33)", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ]))
)[Subtype %in% c("Normal","LumA","LumB","Her2","Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal","LumA","LumB","Her2","Basal"))]
p4b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_boxplot(width = 0.55, outlier.size = 0.45, alpha = 0.88) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal"="#94A3B8","LumA"="#3B82F6","LumB"="#6366F1","Her2"="#F59E0B","Basal"="#E11D48")) +
  labs(title = "B. PAM50 Subtype Expression (TCGA, n=984)", subtitle = "PPARG highest in Normal/LumA; CXCL8 Basal > LumA", x = NULL, y = "log2(RSEM+1)") +
  theme_journal(9.5) + theme(legend.position = "none")

idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp], expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4c <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_tc_p), data = df_tc_p, pval = TRUE, pval.size = 3.2,
  legend.title = "PPARG", legend.labs = c("Low", "High"), palette = c("#E64B35", "#3C5488"),
  title = sprintf("C. TCGA ER+ OS (n=%d)", nrow(df_tc_p)), xlab = "Months", ylab = "OS", ggtheme = theme_journal(9.2)
)$plot + labs(subtitle = "Univ P=0.0135 | Full Adj P=0.032")

idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp], expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4d <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_mb_p), data = df_mb_p, pval = TRUE, pval.size = 3.2,
  legend.title = "PPARG", legend.labs = c("Low", "High"), palette = c("#E64B35", "#3C5488"),
  title = sprintf("D. METABRIC ER+ OS (n=%d)", nrow(df_mb_p)), xlab = "Months", ylab = "OS", ggtheme = theme_journal(9.2)
)$plot + labs(subtitle = "Cont Cox HR=0.894 (P=0.0064)")

pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4e <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_g20_c), data = df_g20_c, pval = TRUE, pval.size = 3.2,
  legend.title = "CXCL8", legend.labs = c("Low", "High"), palette = c("#3C5488", "#E64B35"),
  title = sprintf("E. GSE20685 CXCL8 OS (n=%d)", nrow(df_g20_c)), xlab = "Years", ylab = "OS", ggtheme = theme_journal(9.2)
)$plot + labs(subtitle = "HR=1.251 (P=5.7e-4) | MFS P=0.0016")

idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern], expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := factor(ifelse(expr >= median(expr), "High", "Low"), levels = c("Low", "High"))]
p4f <- survminer::ggsurvplot(
  survfit(Surv(time, status) ~ Group, data = df_mb_c), data = df_mb_c, pval = TRUE, pval.size = 3.2,
  legend.title = "CXCL8", legend.labs = c("Low", "High"), palette = c("#3C5488", "#E64B35"),
  title = sprintf("F. METABRIC ER- OS (n=%d)", nrow(df_mb_c)), xlab = "Months", ylab = "OS", ggtheme = theme_journal(9.2)
)$plot + labs(subtitle = "ER- Log-rank P=0.0158 | Overall Multivar P=7.3e-4")
rm(tc, mb, g20, g42, cl_tc); gc()

# 4G: 5-Cohort Meta & I2 Decomposition
surv_all <- fread(file.path(res_dir, "module4_survival_analysis.csv"))
meta_hr  <- fread(file.path(res_dir, "planA_multicohort_meta_survival.csv"))
het_dt   <- fread(file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))

hr_sub <- rbind(
  meta_hr[analysis == "PPARG_5Cohort_Overall_OS", .(axis = "PPARG Axis", label = sprintf("PPARG 5-Cohort Overall OS (n=%d, I2=%.0f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Pooled Meta")],
  meta_hr[analysis == "PPARG_2Cohort_ERpos_OS",  .(axis = "PPARG Axis", label = sprintf("PPARG 2-Cohort ER+ OS (n=%d, I2=%.0f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Pooled Meta")],
  het_dt[analysis == "All_5_Cohorts_Overall_OS", .(axis = "CXCL8 Heterogeneity Audit", label = sprintf("CXCL8 All 5 Cohorts OS (n=%d, I2=%.1f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Heterogeneous")],
  het_dt[analysis == "LOO_Excluding_TCGA_BRCA",  .(axis = "CXCL8 Heterogeneity Audit", label = sprintf("CXCL8 Microarray 4-Cohort OS (I2=%.0f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Pooled Meta")],
  het_dt[analysis == "GEO_3Cohorts_Affy_OS",     .(axis = "CXCL8 Heterogeneity Audit", label = sprintf("CXCL8 3 Homogeneous Affy GEO (n=%d, I2=%.0f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Zero-Het Meta")],
  meta_hr[analysis == "CXCL8_5Cohort_Progression_RFS_MFS_PFS", .(axis = "CXCL8 Progression", label = sprintf("CXCL8 5-Cohort RFS/MFS/PFS (n=%d, I2=%.1f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Pooled Meta")],
  meta_hr[analysis == "CXCL8_4Ext_RFS_MFS",      .(axis = "CXCL8 Progression", label = sprintf("CXCL8 4-Ext-Cohort RFS/MFS (n=%d, I2=%.1f%%)", total_n, I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Pooled Meta")]
)
hr_sub[, ord := factor(label, levels = rev(label))]
p4g <- ggplot(hr_sub, aes(x = HR, y = ord, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.28, linewidth = 0.78) +
  geom_point(size = 3.2, shape = 18) +
  geom_text(aes(x = HR_UCI95 + 0.015, label = sprintf("HR=%.3f, P=%.2g", HR, pval)),
            hjust = 0, size = 2.75, color = "#0F172A", show.legend = FALSE) +
  scale_color_manual(values = c("Zero-Het Meta" = pal_npg[["red"]], "Pooled Meta" = pal_npg[["blue"]], "Heterogeneous" = "#94A3B8")) +
  coord_cartesian(xlim = c(0.75, 1.56)) +
  labs(title = "G. Meta-Survival & Heterogeneity Decomposition",
       subtitle = "CXCL8 OS I2 drops 78.3% -> 0.0% across Affy GEO (HR=1.248, P=2.11e-6); PPARG ER+ P=3.8e-4",
       x = "Pooled Random-Effects HR (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

# 4H: PPARG Adipose-Confounding & Multivariable Cox Audit
ad_dt <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
ad_sub <- ad_dt[cohort %in% c("TCGA_ERpos", "TCGA_ERpos_LowADIPOQ", "TCGA_ERpos_LowPLIN1",
                              "METABRIC_ERpos", "METABRIC_ERpos_LowADIPOQ", "METABRIC_ERpos_LowPLIN1")]
ad_sub[, coh_grp := ifelse(grepl("TCGA", cohort), "TCGA-BRCA ER+ (RNA-seq, n=770-783)", "METABRIC ER+ (Illumina, n=1,127-1,579)")]
ad_sub[, lbl := fifelse(model == "M1_Univariate", "M1: Univariate PPARG",
                fifelse(model == "M2_Adj_Adipose(ADIPOQ+PLIN1)", "M2: + ADIPOQ + PLIN1",
                fifelse(model == "M3_Adj_Adipose+Stromal", "M3: + Adipose + StromalScore",
                fifelse(model == "M4_Adj_Age+Stage", "M4: + Age + Stage",
                fifelse(model == "M3_Adj_Age+Stage", "M3: + Age + Stage",
                fifelse(model == "M4_Adj_Age+Stage+Grade", "M4: + Age + Stage + Grade",
                fifelse(grepl("M5", model), "M5: Full Multivar (+Adipose+Clin)",
                fifelse(model == "M6_LowAdipose_Subgroup(ADIPOQ<med)", "M6: Low-ADIPOQ Stratum (<median)",
                        "M7: Low-PLIN1 Stratum (<median)"))))))))]
ad_sub[, sig := ifelse(pval < 0.05, "Significant (P < 0.05)", "Attenuated (P >= 0.05)")]
ad_sub[, row_ord := paste0(coh_grp, "___", lbl)]
ad_sub[, row_ord := factor(row_ord, levels = rev(row_ord))]
p4h <- ggplot(ad_sub, aes(x = HR, y = row_ord, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.26, linewidth = 0.72) +
  geom_point(size = 2.6) +
  geom_text(aes(x = HR_UCI95 + 0.015, label = sprintf("%.3f (P=%.3f)", HR, pval)), hjust = 0, size = 2.6, show.legend = FALSE) +
  facet_wrap(~ coh_grp, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_color_manual(values = c("Significant (P < 0.05)" = pal_npg[["blue"]], "Attenuated (P >= 0.05)" = "#94A3B8")) +
  coord_cartesian(xlim = c(0.55, 1.38)) +
  labs(title = "H. Adipose Confounding & Multivariable Cox Audit (PPARG ER+ OS)",
       subtitle = "TCGA ER+ persists across all adipose/stromal/clinical adjustments (M5 P=0.032; Low-ADIPOQ P=0.012); METABRIC attenuates",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.2) + theme(legend.position = "bottom")

# 4I: TCGA ER+ Nomogram ROC/C-index Incremental Gain
roc_dt <- fread(file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))[cohort == "TCGA_ERpos" & model != "M1_PPARG_Alone"]
roc_dt[, Model := factor(model,
                         levels = c("M2_Clinical_Base","M3_Clin_Plus_PPARG","M4_Clin_Adipose_Base","M5_Clin_Adipose_Plus_PPARG"),
                         labels = c("Age+Stage","Age+Stage+PPARG\n(LRT P=0.036)","Clin+Adipose+Stromal","Full+PPARG\n(LRT P=0.028)"))]
roc_m <- melt(roc_dt, id.vars = "Model", measure.vars = c("C_index","AUC_3yr","AUC_5yr","AUC_10yr"),
              variable.name = "Metric", value.name = "Score")
roc_m[, Metric := factor(Metric, levels = c("C_index","AUC_3yr","AUC_5yr","AUC_10yr"),
                         labels = c("C-index","3-yr AUC","5-yr AUC","10-yr AUC"))]
p4i <- ggplot(roc_m, aes(x = Metric, y = Score, fill = Model)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", Score)), position = position_dodge(width = 0.78), vjust = -0.35, size = 2.35) +
  scale_fill_manual(values = c("#94A3B8","#3C5488","#CBD5E1","#00A087")) +
  coord_cartesian(ylim = c(0.60, 0.87)) +
  labs(title = "I. TCGA ER+ Nomogram C-index & timeROC Gain",
       subtitle = "Adding PPARG improves C-index, AIC (922.9->920.0), & 10-yr AUC (0.679->0.699)",
       x = NULL, y = "Discrimination Score", fill = NULL) +
  theme_journal(9.2) + theme(legend.position = "bottom")

# 4J: 5-Yr Calibration Plot from nom_obj
nom_obj <- readRDS(file.path(bulk_dir, "planA_tcga_erpos_nomogram.rds"))
cal_df <- as.data.table(unclass(nom_obj$cal_5yr))
cal_df[, lci := pmax(0, KM.corrected - 1.96 * std.err)]
cal_df[, uci := pmin(1, KM.corrected + 1.96 * std.err)]
p4j <- ggplot(cal_df, aes(x = mean.predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#64748B", linewidth = 0.7) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.015, color = "#3C5488", linewidth = 0.65) +
  geom_line(aes(y = KM, color = "Apparent KM"), linewidth = 0.85, linetype = "dotted") +
  geom_point(aes(y = KM, color = "Apparent KM"), size = 2.4) +
  geom_line(aes(y = KM.corrected, color = "Bootstrap Corrected (B=200)"), linewidth = 1.0) +
  geom_point(aes(y = KM.corrected, color = "Bootstrap Corrected (B=200)"), size = 3.0) +
  scale_color_manual(values = c("Apparent KM" = "#F39B7F", "Bootstrap Corrected (B=200)" = "#3C5488")) +
  coord_fixed(xlim = c(0.42, 1.0), ylim = c(0.42, 1.0)) +
  labs(title = "J. 5-Year OS Nomogram Calibration (TCGA ER+, n=770)",
       subtitle = "Internal bootstrap B=200 validation (Age + Stage + PPARG)",
       x = "Nomogram-Predicted 5-Yr OS", y = "Observed 5-Yr OS (KM)", color = NULL) +
  theme_journal(9.2) + theme(legend.position = "bottom")

fig4_r1 <- (p4a | p4b)
fig4_r2 <- (p4c | p4d | p4e | p4f)
fig4_r3 <- (p4g | p4h)
fig4_r4 <- (p4i | p4j)
fig4_all <- fig4_r1 / fig4_r2 / fig4_r3 / fig4_r4 + plot_layout(heights = c(0.85, 0.90, 1.15, 0.95))
save_pub(fig4_all, "Fig4_Pub_Clinical_Multivar_Nomogram_and_Meta", w = 17.5, h = 21.0)
rm(fig4_all, fig4_r1, fig4_r2, fig4_r3, fig4_r4, p4a, p4b, p4c, p4d, p4e, p4f, p4g, p4h, p4i, p4j, nom_obj); gc()

# ==============================================================================
# FIG 5: 真实乳腺癌单细胞 TISCH2 定位 (5A/B) + HPA IHC & nCPM (5C) + TME (5D) + GSEA (5E)
# ==============================================================================
sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
sc_sum   <- fread(file.path(res_dir, "planA_scrna_celltype_summary.csv"))
# Prioritize GSE176078 if available, otherwise EMTAB8107 cleanly
sc_primary_name <- if ("BRCA_GSE176078" %in% names(sc_cache)) "BRCA_GSE176078" else "BRCA_EMTAB8107"
sc_df   <- copy(sc_cache[[sc_primary_name]])
sc_umap <- sc_df[!is.na(UMAP_1) & !is.na(UMAP_2)]
sc_sum_p <- sc_sum[dataset == sc_primary_name]

pct_lookup <- function(g, ct) {
  v <- sc_sum_p[gene == g & celltype_major == ct, pct_expr]
  if (length(v) == 0) NA_real_ else v[1]
}
overall_pct <- function(g) {
  d <- sc_sum_p[gene == g]
  if (nrow(d) == 0) 0 else round(100 * sum(d$n_expr) / sum(d$n_cells), 2)
}
max_pct <- function(g) {
  d <- sc_sum_p[gene == g]
  if (nrow(d) == 0) 0 else max(d$pct_expr, na.rm = TRUE)
}

# UMAP Color Palette
ct_levels <- sc_umap[, .N, by = Celltype_major][order(-N)]$Celltype_major
ct_pal <- setNames(colorRampPalette(c("#3C5488","#E64B35","#00A087","#4DBBD5","#F39B7F","#8491B4","#D97706","#7C3AED","#EC4899","#14B8A6","#64748B"))(length(ct_levels)), ct_levels)

p5a_ct <- ggplot(sc_umap, aes(x = UMAP_1, y = UMAP_2, color = Celltype_major)) +
  geom_point(size = 0.22, alpha = 0.65) +
  scale_color_manual(values = ct_pal, guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = sprintf("A1. Tumor scRNA (%s, n=%s)", sc_primary_name, format(nrow(sc_df), big.mark=",")),
       subtitle = sprintf("ADIPOQ=%.2f%%; PLIN1<=%.2f%% (UMAP n=%s)", overall_pct("ADIPOQ"), max_pct("PLIN1"), format(nrow(sc_umap), big.mark=",")),
       x = "UMAP 1", y = "UMAP 2", color = "Lineage") +
  theme_journal(8.8) + theme(legend.position = "right", legend.key.size = unit(0.32, "cm"))

plot_umap_gene <- function(dt, gene, ttl, sub_t) {
  d <- copy(dt[, .(UMAP_1, UMAP_2, val = get(gene))])
  d <- d[order(val)] # draw expressing cells on top
  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = val)) +
    geom_point(size = 0.22, alpha = ifelse(d$val > 0, 0.85, 0.25)) +
    scale_color_gradient(low = "#E2E8F0", high = "#DC2626") +
    labs(title = ttl, subtitle = sub_t, x = "UMAP 1", y = "UMAP 2", color = "Expr") +
    theme_journal(8.8) + theme(legend.position = "right", legend.key.height = unit(0.55, "cm"))
}

p5a_cx <- plot_umap_gene(
  sc_umap, "CXCL8", "A2. CXCL8 Tumor scRNA Localization",
  sprintf("Mono/Macro %.1f%% vs Malignant %.1f%%", pct_lookup("CXCL8", "Mono/Macro"), pct_lookup("CXCL8", "Malignant"))
)
p5a_pp <- plot_umap_gene(
  sc_umap, "PPARG", "A3. PPARG Tumor scRNA Localization",
  sprintf("Endo %.1f%%, Fibro %.1f%%, Macro %.1f%%, Malig %.1f%%",
          pct_lookup("PPARG", "Endothelial"), pct_lookup("PPARG", "Fibroblasts"),
          pct_lookup("PPARG", "Mono/Macro"), pct_lookup("PPARG", "Malignant"))
)
p5a_pl <- plot_umap_gene(
  sc_umap, "PLIN1", sprintf("A4. Adipocyte Ctrl: PLIN1 (ADIPOQ %.2f%%)", overall_pct("ADIPOQ")),
  sprintf("Virtually absent (PLIN1 <=%.2f%%, LEP <=%.2f%%)", max_pct("PLIN1"), max_pct("LEP"))
)

# 5B: DotPlot across lineages
dot_genes <- intersect(c("CXCL8","IL1B","FFAR2","PPARG","FABP4","ADIPOQ","PLIN1","LEP","EPCAM","VWF","CD68"), unique(sc_sum$gene))
dot_dt <- sc_sum[dataset == sc_primary_name & gene %in% dot_genes]
dot_dt[, gene := factor(gene, levels = rev(dot_genes))]
ct_ord <- c("Malignant","Epithelial","Endothelial","PVL","Myofibroblasts","Fibroblasts","Mono/Macro","Mast","T_cells","CD8T","CD8Tex","CD4Tconv","Tprolif","B","Plasma","DC","SMC")
dot_dt[, celltype_major := factor(celltype_major, levels = intersect(ct_ord, unique(celltype_major)))]

p5b <- ggplot(dot_dt, aes(x = celltype_major, y = gene, size = pct_expr, color = mean_expr)) +
  geom_point() +
  scale_size(range = c(0.4, 6.0), name = "% Expressing") +
  scale_color_gradient(low = "#E0F2FE", high = "#BE123C", name = "Mean Expr") +
  labs(title = sprintf("B. Single-Cell Lineage Specificity DotPlot (%s)", sc_primary_name),
       subtitle = "CXCL8/IL1B track Mono/Macro; PPARG+FABP4 track Endothelial/PVL while ADIPOQ/PLIN1/LEP stay flat",
       x = NULL, y = NULL) +
  theme_journal(9.0) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

# Subtype panel if GSE176078 present, else HPA nCPM barplot
if ("BRCA_GSE176078" %in% names(sc_cache) && file.exists(file.path(res_dir, "planA_scrna_subtype_comparison.csv"))) {
  st_cmp <- fread(file.path(res_dir, "planA_scrna_subtype_comparison.csv"))
  st_sub <- st_cmp[celltype %in% c("ALL_CELLS", "Endothelial", "Malignant", "Mono/Macro") & gene %in% c("PPARG", "CXCL8")]
  st_sub[, celltype := factor(celltype, levels = c("ALL_CELLS", "Endothelial", "Malignant", "Mono/Macro"),
                              labels = c("All Cells", "Endothelial", "Malignant", "Mono/Macro"))]
  st_sub[, subtype := factor(subtype, levels = c("ER+", "HER2+", "TNBC"))]
  st_sub[, gene := factor(gene, levels = c("PPARG", "CXCL8"))]
  p5c_left <- ggplot(st_sub, aes(x = subtype, y = mean_expr, fill = subtype)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.2f", mean_expr)), vjust = -0.25, size = 2.05) +
    facet_grid(gene ~ celltype, scales = "free_y") +
    scale_fill_manual(values = c("ER+"="#3B82F6","HER2+"="#F59E0B","TNBC"="#E11D48")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = "C1. GSE176078 Subtypes (n=100,064)",
         subtitle = "All/Endo PPARG ER+>TNBC (P<1e-45); Malig CXCL8 TNBC>ER+ (P=6e-63)",
         x = NULL, y = "Mean Norm. Expr") +
    theme_journal(8.4) + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, size = 7.3))
} else {
  sc_hpa <- fread(file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
  sc_keep <- c("Neutrophils","Macrophages","Fibroblasts","Adipocytes","Lymphatic endothelial cells","Vascular endothelial cells","Breast secretory cells","T-cells","B-cells")
  sc_sub <- sc_hpa[gene_symbol %in% c("CXCL8","PPARG") & cell_type %in% sc_keep]
  sc_sub[, log10_ncpm := log10(nCPM + 1)]
  sc_ord <- sc_sub[, .(mx = max(log10_ncpm)), by = cell_type][order(mx)]$cell_type
  sc_sub[, cell_type := factor(cell_type, levels = sc_ord)]
  p5c_left <- ggplot(sc_sub, aes(x = log10_ncpm, y = cell_type, fill = gene_symbol)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.64) +
    scale_fill_manual(values = c("CXCL8" = pal_npg[["red"]], "PPARG" = pal_npg[["blue"]])) +
    labs(title = "C1. HPA Reference Single-Cell nCPM", subtitle = "PPARG Endothelial=234, Adipocyte=2599; CXCL8 Neutrophil=5782", x = "log10(nCPM+1)", y = NULL, fill = NULL) +
    theme_journal(8.8) + theme(legend.position = "bottom")
}

img_panel <- function(jpg_path, ttl, sub_t) {
  img <- jpeg::readJPEG(jpg_path)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + coord_fixed(ratio=1) +
    labs(title = ttl, subtitle = sub_t) + theme_void() +
    theme(aspect.ratio = 1, plot.title = element_text(face="bold", size=8.2, hjust=0.5),
          plot.subtitle = element_text(size=7.2, color="#475569", hjust=0.5),
          panel.border = element_rect(color="#CBD5E1", fill=NA, linewidth=0.55))
}
p5c_ihc <- (
  img_panel(file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"), "C2. Normal PPARG", "Adipocyte Low|Duct ND") |
  img_panel(file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"), "C3. BRCA PPARG IHC", "4 Med/6 Low/2 ND") |
  img_panel(file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"), "C4. BRCA CXCL8 IHC", "12/12 ND (Secreted)")
)

# 5D & 5E: TME & GSEA
cor5 <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
tme_top <- cor5[((gene_symbol=="CXCL8" & stratum=="ER_neg") | (gene_symbol=="PPARG" & stratum=="ER_pos")) &
                  feature_type %in% c("TME_Cell_Score","Checkpoint","Genomic_Feature") & FDR < 0.05][
                    order(gene_symbol, -abs(rho))][, .SD[1:min(10, .N)], by = gene_symbol]
tme_top[, feat_clean := sub("^(quanTIseq|MCPcounter|ESTIMATE)__(.*)$", "\\2 [\\1]", feature)]
tme_top[, facet_lbl := ifelse(gene_symbol=="CXCL8", "D1. CXCL8 in ER- BRCA (n=237)", "D2. PPARG in ER+ BRCA (n=783)")]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feat_clean)]
tme_top <- tme_top[order(facet_lbl, rho)]; tme_top[, feat_ord := factor(feat_ord, levels = unique(feat_ord))]
p5d <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_segment(aes(x = 0, xend = rho, y = feat_ord, yend = feat_ord), linewidth = 0.75) +
  geom_point(size = 2.5) + facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("Positive" = pal_npg[["red"]], "Negative" = pal_npg[["blue"]])) +
  labs(title = "D. Immune Microenvironment & Genomic Correlations (FDR<0.05)", x = "Spearman rho", y = NULL, color = NULL) +
  theme_journal(9.0) + theme(legend.position = "bottom")

gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene=="CXCL8" & stratum=="ER_neg") | (target_gene=="PPARG" & stratum=="ER_pos")) & p.adjust < 0.05][
  order(target_gene, -abs(NES))][, .SD[1:min(9, .N)], by = target_gene]
gs_top[, clean_id := gsub("^HALLMARK_", "", ID)]
gs_top[, facet_lbl := ifelse(target_gene=="CXCL8", "E1. CXCL8 High vs Low (ER- BRCA)", "E2. PPARG High vs Low (ER+ BRCA)")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]; gs_top[, id_ord := factor(id_ord, levels = id_ord)]
p5e <- ggplot(gs_top, aes(x = NES, y = id_ord, fill = -log10(p.adjust))) +
  geom_col(width = 0.72) + facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_fill_gradient(low = "#93C5FD", high = "#1E3A8A") +
  labs(title = "E. Single-Gene Hallmark GSEA (FDR<0.05)", x = "Normalized Enrichment Score (NES)", y = NULL, fill = "-log10(FDR)") +
  theme_journal(9.0) + theme(legend.position = "bottom")

fig5_r1 <- (p5a_ct | p5a_cx | p5a_pp | p5a_pl)
fig5_r2 <- (p5b | p5c_left | p5c_ihc) + plot_layout(widths = c(1.05, 1.05, 1.05))
fig5_r3 <- (p5d | p5e)
fig5 <- fig5_r1 / fig5_r2 / fig5_r3 + plot_layout(heights = c(0.95, 1.0, 1.05))
save_pub(fig5, "Fig5_Pub_Tumor_scRNA_HPA_TME_GSEA", w = 18.0, h = 15.2)
rm(fig5, fig5_r1, fig5_r2, fig5_r3, sc_cache, sc_df, sc_umap, p5a_ct, p5a_cx, p5a_pp, p5a_pl, p5b, p5c_left, p5c_ihc, p5d, p5e); gc()

# ==============================================================================
# FIG 6: 真实临床治疗队列 (6A内分泌 / 6B新辅助化疗DRFS / 6C I-SPY2交互) + GDSC2 细胞系 (6D/6E)
# ==============================================================================
th_dt <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))
th_pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

# 6A: Endocrine Therapy Specificity Forest (METABRIC HT YES vs NO + GSE9195 + GSE17705)
endo_dt <- th_dt[gene == "PPARG" & cohort %in% c("METABRIC", "GSE9195", "GSE17705") & stratum == "ER_pos"]
endo_dt[, lbl := sprintf("%s | %s %s (n=%d)", cohort, gsub("_", " ", treatment_arm), endpoint, n_samples)]
endo_dt[, tag := ifelse(pval < 0.05 | km_or_wilcox_pval < 0.05, "P < 0.05", "Null")]
endo_dt[, lbl := factor(lbl, levels = rev(lbl))]
p6a <- ggplot(endo_dt, aes(x = estimate, y = lbl, color = tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = lci95, xmax = uci95), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(size = 2.8) +
  geom_text(aes(x = pmin(1.52, uci95 * 1.04), label = sprintf("HR=%.3f, P=%.3g", estimate, pval)), hjust = 0, size = 2.55, show.legend = FALSE) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1.0, 1.5)) +
  scale_color_manual(values = c("P < 0.05" = pal_npg[["blue"]], "Null" = "#94A3B8")) +
  coord_cartesian(xlim = c(0.075, 3.6)) +
  labs(title = "A. ER+ Endocrine Therapy Outcomes (PPARG)",
       subtitle = "METABRIC HT=YES: OS P=0.0065, RFS P=0.029 (HT=NO null); GSE9195 DMFS KM P=0.042",
       x = "Hazard Ratio (95% CI, log scale)", y = NULL, color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

# 6B: NACT pCR (GSE25066 + GSE32646) + GSE25066 DRFS KM curve
nact_pcr_dt <- th_dt[gene == "CXCL8" & cohort %in% c("GSE25066", "GSE32646") & endpoint == "pCR_vs_RD" & stratum %in% c("Overall", "ER_neg", "TNBC")]
nact_pcr_dt[, lbl := sprintf("%s %s (n=%d, %d pCR)", cohort, stratum, n_samples, n_events)]
nact_pcr_dt[, tag := ifelse(pval < 0.05, "P < 0.05 (Unadj)", "ER-Stratified")]
nact_pcr_dt[, lbl := factor(lbl, levels = rev(lbl))]
p6b_pcr <- ggplot(nact_pcr_dt, aes(x = estimate, y = lbl, color = tag)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = lci95, xmax = uci95), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(size = 2.6) +
  geom_text(aes(x = uci95 + 0.02, label = sprintf("OR=%.2f, AUC=%.2f, P=%.3g", estimate, auc, pval)), hjust = 0, size = 2.45, show.legend = FALSE) +
  scale_color_manual(values = c("P < 0.05 (Unadj)" = pal_npg[["red"]], "ER-Stratified" = "#64748B")) +
  coord_cartesian(xlim = c(0.55, 2.95)) +
  labs(title = "B. NACT pCR: GSE25066 & GSE32646 (CXCL8)",
       subtitle = "Overall NACT pCR P=0.0083 / 0.0209 (ER-mediated; ER-adj P>0.21)",
       x = "Odds Ratio for pCR (95% CI)", y = NULL, color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

g25_km <- copy(th_pc$gse25066)[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs)]
g25_km[, Group := factor(ifelse(CXCL8 >= median(CXCL8), "High CXCL8", "Low CXCL8"), levels = c("Low CXCL8", "High CXCL8"))]
p6b <- survminer::ggsurvplot(
  survfit(Surv(t_drfs, e_drfs) ~ Group, data = g25_km), data = g25_km, pval = TRUE, pval.size = 3.1,
  legend.title = NULL, palette = c("#3C5488", "#E64B35"),
  title = sprintf("C. GSE25066 NACT DRFS (n=%d)", nrow(g25_km)),
  xlab = "Years", ylab = "Distant Relapse-Free Survival", ggtheme = theme_journal(8.8)
)$plot + labs(subtitle = "Univ HR=1.187 (P=4.7e-5) | ER+pCR Adj HR=1.146 (P=0.003)")

# 6C: I-SPY2 (GSE173839) Durvalumab/Olaparib vs Control Arm Treatment-Interaction!
g173_df <- copy(th_pc$gse173839)[!is.na(pCR)]
g173_df[, Resp := factor(ifelse(pCR == 1, "pCR", "RD"), levels = c("RD", "pCR"))]
g173_df[, Arm_lbl := ifelse(Arm == "durvalumab/olaparib", "Durvalumab + Olaparib\n(n=71; 29 pCR)", "Control Paclitaxel\n(n=29; 7 pCR)")]
p6c1 <- ggplot(g173_df, aes(x = Resp, y = CXCL8, fill = Resp)) +
  geom_boxplot(width = 0.48, alpha = 0.85, outlier.size = 0.8) +
  geom_jitter(width = 0.14, size = 1.1, alpha = 0.55) +
  facet_wrap(~ Arm_lbl) +
  scale_fill_manual(values = c("RD" = "#94A3B8", "pCR" = "#E64B35")) +
  labs(title = "D1. I-SPY2 CXCL8 x Immunotherapy (P_int=0.0056)",
       subtitle = "Durva/Ola: OR=1.73, P=0.0029, AUC=0.712 | Control: OR=0.65", x = NULL, y = "CXCL8 Expression") +
  theme_journal(8.8) + theme(legend.position = "none")

g173_tn <- g173_df[HR_neg == TRUE]
g173_tn[, Arm_lbl := ifelse(Arm == "durvalumab/olaparib", "TNBC Durva+Ola\n(n=21; 9 pCR)", "TNBC Control\n(n=17; 6 pCR)")]
p6c2 <- ggplot(g173_tn, aes(x = Resp, y = PPARG, fill = Resp)) +
  geom_boxplot(width = 0.48, alpha = 0.85, outlier.size = 0.8) +
  geom_jitter(width = 0.14, size = 1.1, alpha = 0.55) +
  facet_wrap(~ Arm_lbl) +
  scale_fill_manual(values = c("RD" = "#94A3B8", "pCR" = "#3C5488")) +
  labs(title = "D2. I-SPY2 PPARG x Arm in TNBC (P_int=0.049)",
       subtitle = "Low PPARG predicts Durva/Ola pCR (P=0.020, AUC=0.861)", x = NULL, y = "PPARG Expression") +
  theme_journal(8.8) + theme(legend.position = "none")

# 6D & 6E: GDSC2 Pharmacogenomics (n=47 BRCA cell lines, Subtype-Adjusted)
sub_d <- fread(file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
d_dedup <- sub_d[order(gene_symbol, all_spearman_pval)][, .SD[1], by = .(gene_symbol, DRUG_NAME)]
d_top <- d_dedup[all_spearman_pval < 0.05][order(gene_symbol, all_spearman_pval)][, .SD[1:min(9, .N)], by = gene_symbol]
d_top[, dir := fifelse(subtype_adjusted_pval < 0.05,
               fifelse(all_spearman_rho < 0, "Sensitive + Subtype-Adj P<0.05", "Resistant + Subtype-Adj P<0.05"),
               fifelse(all_spearman_rho < 0, "Sensitive (Overall P<0.05)", "Resistant (Overall P<0.05)"))]
d_top[, drug_ord := paste0(gene_symbol, "___", DRUG_NAME)]
d_top <- d_top[order(gene_symbol, -all_spearman_pval)]; d_top[, drug_ord := factor(drug_ord, levels = drug_ord)]
p6d <- ggplot(d_top, aes(x = all_spearman_rho, y = drug_ord, color = dir, size = -log10(all_spearman_pval))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B") +
  geom_point() + facet_wrap(~ gene_symbol, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c(
    "Sensitive + Subtype-Adj P<0.05" = "#059669", "Sensitive (Overall P<0.05)" = "#6EE7B7",
    "Resistant + Subtype-Adj P<0.05" = "#DC2626", "Resistant (Overall P<0.05)" = "#FCA5A5"
  )) +
  labs(title = "E. GDSC2 Drug Sensitivity (n=47 BRCA Lines)",
       subtitle = "Fulvestrant/Vorinostat/POMHEX survive subtype OLS adjustment (P<0.05)",
       x = "Spearman rho (Expr vs LN_IC50)", y = NULL, color = NULL, size = "-log10(P)") +
  theme_journal(8.8) + theme(legend.position = "bottom")

gb_mat <- readRDS(file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))
p6e1 <- ggplot(gb_mat[DRUG_NAME == "Dasatinib"], aes(x = CXCL8, y = LN_IC50, color = Subtype_Group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75) +
  geom_smooth(data = gb_mat[DRUG_NAME == "Dasatinib" & Subtype_Group == "TNBC"], method = "lm", se = FALSE, color = "#E11D48", linewidth = 0.8) +
  scale_color_manual(values = c("ER_pos"="#3B82F6","HER2_pos"="#F59E0B","TNBC"="#E11D48")) +
  labs(title = "F1. CXCL8 vs Dasatinib (rho=-0.466, P=0.0011)", subtitle = "Within TNBC (n=21): rho=-0.417, P=0.060",
       x = "CXCL8 (log2 TPM+1)", y = "LN_IC50", color = "Subtype") +
  theme_journal(8.8) + theme(legend.position = "bottom")

p6e2 <- ggplot(gb_mat[DRUG_ID == 1200], aes(x = PPARG, y = LN_IC50, color = Subtype_Group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#0F172A", linetype = "dashed", linewidth = 0.75) +
  geom_smooth(data = gb_mat[DRUG_ID == 1200 & Subtype_Group == "ER_pos"], method = "lm", se = FALSE, color = "#3B82F6", linewidth = 0.8) +
  scale_color_manual(values = c("ER_pos"="#3B82F6","HER2_pos"="#F59E0B","TNBC"="#E11D48")) +
  labs(title = "F2. PPARG vs Fulvestrant (adj P=0.0096)", subtitle = "All (n=46): rho=+0.375; ER+ (n=14): rho=+0.503",
       x = "PPARG (log2 TPM+1)", y = "LN_IC50", color = "Subtype") +
  theme_journal(8.8) + theme(legend.position = "bottom")

fig6 <- (p6a | p6b_pcr | p6b) / (p6c1 | p6c2) / (p6d | p6e1 | p6e2) + plot_layout(heights = c(0.95, 0.90, 1.05))
save_pub(fig6, "Fig6_Pub_Clinical_Therapy_and_GDSC2", w = 18.0, h = 15.2)
rm(fig6, p6a, p6b_pcr, p6b, p6c1, p6c2, p6d, p6e1, p6e2, gb_mat, th_pc); gc()

# ==============================================================================
# FIG S2: rms::nomogram 图表 + CXCL8 10-层异质性分解 + 菌群跨队列 3-Cohort 森林图 + GSE9195 KM
# ==============================================================================
# S2A: Render classic rms::nomogram cleanly to temporary PNG and embed
tc_ad <- readRDS(file.path(bulk_dir, "tcga_planA_adipose_clin.rds"))
d_nom <- as.data.frame(tc_ad[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event) & !is.na(AGE) & !is.na(STAGE) & !is.na(PPARG),
                             .(OS_MONTHS, OS_event, AGE, STAGE = factor(STAGE), PPARG = round(PPARG, 2))])
dd_nom <- datadist(d_nom)
assign("dd_nom", dd_nom, envir = .GlobalEnv)
options(datadist = "dd_nom")
fit_nom <- cph(Surv(OS_MONTHS, OS_event) ~ AGE + STAGE + PPARG, data = d_nom, x = TRUE, y = TRUE, surv = TRUE)
sv_nom  <- Survival(fit_nom)
nom_plot <- nomogram(fit_nom,
                     fun = list(function(x) sv_nom(36, x), function(x) sv_nom(60, x), function(x) sv_nom(120, x)),
                     funlabel = c("3-Year OS", "5-Year OS", "10-Year OS"), lp = FALSE)
tmp_nom_png <- file.path(pub_dir, "_tmp_nomogram.png")
png(tmp_nom_png, width = 2200, height = 1150, res = 220, bg = "white")
par(mar = c(2.5, 2.0, 2.5, 2.0))
plot(nom_plot, xfrac = 0.25, cex.axis = 0.78, cex.var = 0.85, lmgp = 0.25)
title("A. TCGA-BRCA ER+ (n=770) Clinical-PPARG Nomogram (LRT P=0.036)", cex.main = 1.05)
dev.off()
nom_img <- png::readPNG(tmp_nom_png)
ps2a <- ggplot() + annotation_custom(grid::rasterGrob(nom_img, interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void() + theme(panel.border = element_rect(color="#CBD5E1", fill=NA, linewidth=0.5))
unlink(tmp_nom_png)

# S2B: Complete CXCL8 OS Meta Heterogeneity Audit (All 10 strata)
het_dt[, lbl := factor(analysis, levels = rev(analysis))]
ps2b <- ggplot(het_dt, aes(x = RE_HR, y = lbl, color = I2_pct == 0)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = RE_LCI, xmax = RE_UCI), orientation = "y", width = 0.26, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(x = RE_UCI + 0.012, label = sprintf("RE=%.3f, I2=%.1f%%, P=%.2g", RE_HR, I2_pct, RE_pval)),
            hjust = 0, size = 2.45, color = "#0F172A", show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = pal_npg[["red"]], "FALSE" = pal_npg[["blue"]]),
                     labels = c("TRUE" = "I2 = 0% Homogeneous", "FALSE" = "I2 > 0%")) +
  coord_cartesian(xlim = c(0.88, 1.62)) +
  labs(title = "B. CXCL8 OS Meta Leave-One-Out & Platform Stratification",
       subtitle = "Pooling homogeneous Affymetrix GEO eliminates heterogeneity (I2=78.3% -> 0.0%, Q P=0.54)",
       x = "Random-Effects HR (95% CI)", y = NULL, color = NULL) +
  theme_journal(8.8) + theme(legend.position = "bottom")

# S2C: Cross-Cohort Microbiome Replication Forest (Core 6 Taxa across MiBioGen + DMP + FR02 + Meta)
cc_meta <- fread(file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))
cc_rep  <- fread(file.path(res_dir, "module3_dmp_fr02_core6_replication.csv"))
mib_core <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))[gwas_name != "Marvinbryantia"]

# Focus on the 4 core taxa-outcome pairs highlighted in Plan A + Meta
core_pairs <- rbind(
  data.table(pair="Dorea -> ER-", coh="MiBioGen (Disc.)", OR=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR][1], lci=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=mib_core[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Dorea -> ER-", coh="DMP_2022 (Repl.)", OR=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Dorea" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Dorea -> ER-", coh="Pooled FE Meta",   OR=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_OR][1], lci=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_LCI][1], uci=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_UCI][1], p=cc_meta[core_taxon=="Dorea" & outcome_label=="BRCA_ER_neg", FE_pval][1]),
  data.table(pair="Desulfovibrio -> ER-", coh="MiBioGen (Disc.)", OR=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR][1], lci=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=mib_core[gwas_name=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio -> ER-", coh="DMP_2022 (D. piger)", OR=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="DMP_2022" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio -> ER-", coh="FR02_2022 (D. piger)", OR=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR][1], lci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR_LCI95][1], uci=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", OR_UCI95][1], p=cc_rep[gwas_name=="Desulfovibrio piger" & cohort=="FR02_2022" & outcome_label=="BRCA_ER_neg", pval][1]),
  data.table(pair="Desulfovibrio -> ER-", coh="Pooled FE Meta (k=3)", OR=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_OR][1], lci=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_LCI][1], uci=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_UCI][1], p=cc_meta[core_taxon=="Desulfovibrio" & outcome_label=="BRCA_ER_neg", FE_pval][1]),
  data.table(pair="Adlercreutzia -> ER+", coh="MiBioGen (Disc.)", OR=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR][1], lci=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_LCI95][1], uci=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_UCI95][1], p=mib_core[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", pval][1]),
  data.table(pair="Adlercreutzia -> ER+", coh="DMP_2022 (Repl.)", OR=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR][1], lci=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_LCI95][1], uci=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", OR_UCI95][1], p=cc_rep[gwas_name=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", pval][1]),
  data.table(pair="Adlercreutzia -> ER+", coh="Pooled FE Meta", OR=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_OR][1], lci=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_LCI][1], uci=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_UCI][1], p=cc_meta[core_taxon=="Adlercreutzia" & outcome_label=="BRCA_ER_pos", FE_pval][1])
)
core_pairs[, is_meta := grepl("Meta", coh)]
core_pairs[, row_ord := factor(paste(pair, coh, sep=" | "), levels = rev(paste(pair, coh, sep=" | ")))]
ps2c <- ggplot(core_pairs, aes(x = OR, y = row_ord, color = is_meta)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B") +
  geom_errorbar(aes(xmin = lci, xmax = uci), orientation = "y", width = 0.25, linewidth = 0.72) +
  geom_point(aes(shape = is_meta), size = 2.6) +
  geom_text(aes(x = uci + 0.015, label = sprintf("OR=%.3f, P=%.3g", OR, p)), hjust = 0, size = 2.45, color = "#0F172A", show.legend = FALSE) +
  facet_wrap(~ pair, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
  scale_color_manual(values = c("FALSE" = pal_npg[["blue"]], "TRUE" = pal_npg[["red"]]), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  coord_cartesian(xlim = c(0.65, 1.32)) +
  labs(title = "C. Cross-Cohort Gut Microbiome Replication (MiBioGen + DMP + FR02)",
       subtitle = "Consistent directional protection across independent European microbiome cohorts",
       x = "Odds Ratio (95% CI)", y = NULL) +
  theme_journal(8.8)

# S2D: GSE9195 Tamoxifen DMFS KM curve
th_pc2 <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))
g91_km <- copy(th_pc2$gse9195)[!is.na(t_dmfs) & t_dmfs > 0 & !is.na(e_dmfs)]
g91_km[, Group := factor(ifelse(PPARG >= median(PPARG), "High PPARG", "Low PPARG"), levels = c("Low PPARG", "High PPARG"))]
ps2d <- survminer::ggsurvplot(
  survfit(Surv(t_dmfs, e_dmfs) ~ Group, data = g91_km), data = g91_km, pval = TRUE, pval.size = 3.2,
  legend.title = NULL, palette = c("#E64B35", "#3C5488"),
  title = sprintf("D. GSE9195 Tamoxifen ER+ DMFS (n=%d)", nrow(g91_km)),
  xlab = "Years", ylab = "Distant Metastasis-Free Survival", ggtheme = theme_journal(8.8)
)$plot + labs(subtitle = "Log-rank P=0.0421 | Cont HR=0.317 (P=0.083)")

figs2 <- (ps2a | ps2b) / (ps2c | ps2d) + plot_layout(heights = c(1.0, 1.1))
save_pub(figs2, "FigS2_Pub_NomogramRuler_HeterogeneityLOO_DMPFR02", w = 16.8, h = 13.5)

cat("\nALL FINAL PLAN A FIGURES GENERATED SUCCESSFULLY!\n")
