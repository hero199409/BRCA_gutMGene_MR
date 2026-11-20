# 23_generate_frontiers_part1_figs1to5.R
# Nature / Cell / Frontiers 顶刊美学重塑版 Part 1 (Fig 1 - Fig 5)
# 特点：零副标题小字、零代码下划线/ASCII符号、开放式L型主题、圆角卡片流程图、平滑贝塞尔双轴网络、语义严格染色

source(if (file.exists("scripts/23_0_plot_helpers.R")) "scripts/23_0_plot_helpers.R" else "23_0_plot_helpers.R")

# ==============================================================================
# FIG 1: 研究设计 & 多阶段质控漏斗 (1A 圆角卡片) + 上游菌群/代谢物 MR (1B/1C)
# ==============================================================================
cards_1a <- data.table(
  id     = c("s1", "qc1", "s2", "qc2", "s3", "qc3", "qc4", "ax_erneg", "ax_erpos"),
  xmin   = c( 2,  2, 35.5, 35.5, 69, 69, 69,  4, 53),
  xmax   = c(31, 31, 64.5, 64.5, 98, 98, 98, 47, 96),
  ymin   = c(50, 34, 50,   34,   50, 34, 21,  2,  2),
  ymax   = c(61, 45, 61,   45,   61, 45, 31, 13, 13),
  fill   = c("#EFF6FF","#F8FAFC","#ECFDF5","#FFFBEB","#F5F3FF","#F8FAFC","#FFF7ED","#FFF1F2","#EFF6FF"),
  bcol   = c("#BFDBFE","#CBD5E1","#A7F3D0","#FDE68A","#DDD6FE","#CBD5E1","#FED7AA","#FECDD3","#BFDBFE"),
  badge  = c("STAGE 1","QC 1","STAGE 2","QC 2","STAGE 3","QC 3","QC 4","ER\u2212 AXIS","ER+ AXIS"),
  bg_col = c("#3C5488","#64748B","#00A087","#D97706","#6D28D9","#64748B","#EA580C","#E64B35","#3C5488"),
  title  = c(
    "Gut Microbiome Two-Sample MR",
    "Reverse MR & External Meta",
    "Plasma Metabolome MR",
    "Directionality & Mediation",
    "Host Gene eQTL Discovery",
    "Strict \u00b11 Mb cis-IV Filter",
    "Colocalization & Plasma pQTLs",
    "CXCL8 Inflammatory Risk Axis",
    "PPARG Nuclear Receptor Axis"
  ),
  body   = c(
    "64 aligned taxa (192 tests) \u2192 6 core taxa",
    "Excl. Marvinbryantia | DMP + FR02 Meta",
    "47 metabolites (141 tests) \u2192 5 prioritized",
    "Reverse MR P > 0.05 | Step-1 taxa \u2192 met",
    "41 eQTLGen genes (123 tests) \u2192 8 cis-hits",
    "Excl. 8 trans-only \u2192 29 cis-genes",
    "Excl. CXCR2 (PP.H3 = 1.0) | 5 pQTLs",
    "Microbes / Acetate \u2192 CXCL8 (+IL1B / FFAR2)\nMacrophage scRNA | DRFS | I-SPY2 Pint = 0.0056",
    "Adlercreutzia / Linoleic acid \u2192 PPARG\nEndothelial scRNA | Adipose-adj. OS | Endocrine"
  ),
  txt_col = c("#1E3A8A","#334155","#065F46","#92400E","#5B21B6","#334155","#9A3412","#9F1239","#1E3A8A")
)

# 垂直流线与底部直线分支连接（箭头均精准停在胶囊徽章上沿之外）
flows_v <- data.table(
  x1 = c(16.5, 50.0, 83.5, 83.5),
  y1 = c(50.0, 50.0, 50.0, 34.0),
  x2 = c(16.5, 50.0, 83.5, 83.5),
  y2 = c(46.5, 46.5, 46.5, 32.5)
)
flows_diag <- data.table(
  x1 = c(16.5, 50.0, 50.0, 83.5, 83.5),
  y1 = c(34.0, 34.0, 34.0, 21.0, 21.0),
  x2 = c(22.5, 25.5, 72.5, 28.5, 76.5),
  y2 = c(14.5, 14.5, 14.5, 14.5, 14.5)
)

p1a <- ggplot() +
  geom_segment(data = flows_diag, aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "#94A3B8", linewidth = 0.80, alpha = 0.85,
               arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  geom_segment(data = flows_v, aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "#64748B", linewidth = 0.85,
               arrow = arrow(length = unit(0.17, "cm"), type = "closed"))

for (i in seq_len(nrow(cards_1a))) {
  r <- cards_1a[i]
  p1a <- add_round_card(p1a, r$xmin, r$xmax, r$ymin, r$ymax, fill = r$fill, border_col = r$bcol, lwd = 0.8, r_pt = 8)
  bw <- nchar(r$badge) * 0.72 + 2.8
  bx1 <- (r$xmin + r$xmax) / 2 - bw / 2
  bx2 <- (r$xmin + r$xmax) / 2 + bw / 2
  p1a <- add_round_card(p1a, bx1, bx2, r$ymax - 1.4, r$ymax + 1.3, fill = r$bg_col, border_col = NA, r_pt = 6)
}

cards_1a[, `:=`(cx = (xmin + xmax) / 2,
                ty = ymax - 3.3,
                by = (ymin + ymax - 2.2) / 2,
                badge_y = ymax - 0.05)]

p1a <- p1a +
  geom_text(data = cards_1a, aes(x = cx, y = badge_y, label = badge), color = "white", size = 2.25, fontface = "bold") +
  geom_text(data = cards_1a, aes(x = cx, y = ty, label = title, color = txt_col), size = 3.05, fontface = "bold") +
  geom_text(data = cards_1a, aes(x = cx, y = by, label = body), color = "#334155", size = 2.65, lineheight = 1.22) +
  scale_color_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0.5, 64), expand = FALSE) +
  labs(title = "A. Study Design & Multi-Stage QC Funnel") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11.2, color = "#0F172A", margin = margin(b = 6)))

# 1B: 菌群 MR（清理 genus/family 与下划线）
mib_sig <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))[gwas_name != "Marvinbryantia"]
mib_sig[, rank_tag := ifelse(gwas_rank == "family", " (Fam.)", "")]
mib_sig[, label := sprintf("%s%s \u2192 %s (%d SNPs)", gwas_name, rank_tag, clean_subtype_str(outcome_label), nsnp)]
mib_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", fifelse(FDR < 0.10, "FDR < 0.10", "Nominal P < 0.05"))]
mib_sig[, domain_flag := fifelse(nsnp <= 2, "Warning (2 SNPs)", "Primary")]
mib_sig[, anno_txt := sprintf("OR = %.2f [%.2f\u2013%.2f], P = %.2g", OR, OR_LCI95, OR_UCI95, pval)]
mib_sig <- mib_sig[order(pval)]
mib_sig[, label := factor(label, levels = rev(label))]

p1b <- plot_forest_panel(
  mib_sig, "B. Gut Microbiota \u2192 BCAC",
  c(0.62, max(mib_sig$OR_UCI95) * 1.54),
  c("FDR < 0.05" = pal_npg[["red"]], "FDR < 0.10" = pal_npg[["amber"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

met_sig <- fread(file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv"))
met_sig <- met_sig[id.exposure %in% c("met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-446")]
met_sig[, trait_short := fifelse(grepl("hydroxyphen", trait_name), "3-(4-hydroxyphenyl)lactate",
                          fifelse(grepl("phenylpropionate", trait_name), "3-phenylpropionate",
                          fifelse(grepl("linoleic", trait_name), "18:2 Linoleic acid", trait_name)))]
met_sig[, label := sprintf("%s \u2192 %s (%d SNPs)", trait_short, clean_subtype_str(outcome_label), nsnp)]
met_sig[, sig_tag := fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05")]
met_sig[, domain_flag := "Primary"]
met_sig[, anno_txt := sprintf("OR = %.2f [%.2f\u2013%.2f], P = %.2g", OR, OR_LCI95, OR_UCI95, pval)]
met_sig <- met_sig[order(pval)]
met_sig[, label := factor(label, levels = rev(label))]

p1c <- plot_forest_panel(
  met_sig, "C. Plasma Metabolites \u2192 BCAC",
  c(0.45, max(met_sig$OR_UCI95) * 1.56),
  c("FDR < 0.05" = pal_npg[["red"]], "Nominal P < 0.05" = pal_npg[["blue"]])
)

fig1_all <- p1a / (p1b | p1c) + plot_layout(heights = c(1.12, 1.0))
save_pub(fig1_all, "Fig1_Design_QC_and_Upstream_MR", w = 17.5, h = 12.6)
rm(fig1_all, p1a, p1b, p1c, cards_1a, flows_v, flows_diag); gc()

# ==============================================================================
# FIG 2: 纯白画布 + 平滑 S 形贝塞尔曲线 + 圆角胶囊网络图
# ==============================================================================
edges_dt <- fread(file.path(res_dir, "module3_candidate_gutmgene_edges.csv"))
cis_pos  <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_g    <- unique(cis_pos$gene_symbol)

sub_e <- edges_dt[gene_symbol %in% pos_g & !is.na(metabolite)]
map_microbe <- function(m) {
  if (grepl("Parabacteroides", m, ignore.case = TRUE)) return("Parabacteroides (Prot.)")
  if (grepl("Dorea", m, ignore.case = TRUE)) return("Dorea (ER\u2212 Meta)")
  if (grepl("Ruminococc", m, ignore.case = TRUE)) return("Ruminococcaceae (Risk)")
  if (grepl("Blautia", m, ignore.case = TRUE)) return("Blautia (ER\u2212, 2 SNPs)")
  if (grepl("MR_Positive_Metabolite", m)) return("Blood Acetate (FDR < 0.05)")
  return(NA_character_)
}
sub_e[, m_grp := sapply(microbe, map_microbe)]
mm_links <- unique(rbind(
  unique(sub_e[!is.na(m_grp), .(from = m_grp, to = metabolite)]),
  data.table(from = "Adlercreutzia (ER+ Meta)", to = "Equol"),
  data.table(from = "Blood Acetate (FDR < 0.05)", to = "Acetate"),
  data.table(from = "Linoleic Acid (18:2, ER\u2212 MR)", to = "10-keto-12Z-octadecenoic acid")
))
mg_links <- unique(sub_e[, .(from = metabolite, to = gene_symbol, alteration)])
mg_links <- mg_links[, .(alteration = ifelse(length(unique(alteration)) > 1, "mixed", alteration[1])), by = .(from, to)]
mg_links <- unique(rbind(mg_links, data.table(from = "10-keto-12Z-octadecenoic acid", to = "PPARG", alteration = "activation")))

left_nodes  <- c("Parabacteroides (Prot.)", "Dorea (ER\u2212 Meta)", "Blautia (ER\u2212, 2 SNPs)",
                 "Blood Acetate (FDR < 0.05)", "Ruminococcaceae (Risk)",
                 "Linoleic Acid (18:2, ER\u2212 MR)", "Adlercreutzia (ER+ Meta)")
mid_nodes   <- c("Butyrate", "Acetate", "Succinate", "3-Indolepropionic acid",
                 "Trimethylamine oxide", "10-keto-12Z-octadecenoic acid", "Equol")
right_nodes <- c("CXCL8 (ER\u2212 eQTL + 2\u00d7 pQTL)", "IL1B (Coloc PP.H4 = 0.507)", "FFAR2 (GTEx Breast)",
                 "TNF (Secondary)", "HDAC11 (Secondary)", "IL12A (Secondary)",
                 "PPARG (ER+ eQTL + OS)", "CXCR2 (Excl.: PP.H3 = 1.0)")

rename_g <- c(
  "CXCL8"  = "CXCL8 (ER\u2212 eQTL + 2\u00d7 pQTL)",
  "IL1B"   = "IL1B (Coloc PP.H4 = 0.507)",
  "FFAR2"  = "FFAR2 (GTEx Breast)",
  "TNF"    = "TNF (Secondary)",
  "HDAC11" = "HDAC11 (Secondary)",
  "IL12A"  = "IL12A (Secondary)",
  "PPARG"  = "PPARG (ER+ eQTL + OS)",
  "CXCR2"  = "CXCR2 (Excl.: PP.H3 = 1.0)"
)
mg_links[, to := rename_g[to]]
mm_links <- mm_links[from %in% left_nodes & to %in% mid_nodes]
mg_links <- mg_links[from %in% mid_nodes & to %in% right_nodes]

node_pos <- rbindlist(list(
  data.table(node = left_nodes,  col = 1.0, hw = 0.28, y = seq(9.5, 1.5, length.out = length(left_nodes))),
  data.table(node = mid_nodes,   col = 2.0, hw = 0.25, y = seq(9.5, 1.5, length.out = length(mid_nodes))),
  data.table(node = right_nodes, col = 3.0, hw = 0.28, y = seq(9.5, 1.5, length.out = length(right_nodes)))
))
e1 <- merge(merge(mm_links, node_pos[, .(from = node, x1 = col + hw, y1 = y)], by = "from"),
            node_pos[, .(to = node, x2 = col - hw, y2 = y)], by = "to")
e1[, `:=`(edge_type = "Produces", is_core = grepl("Dorea|Parabacteroides|Acetate|Linoleic|Adlercreutzia", from))]
e2 <- merge(merge(mg_links, node_pos[, .(from = node, x1 = col + hw, y1 = y)], by = "from"),
            node_pos[, .(to = node, x2 = col - hw, y2 = y)], by = "to")
e2[, edge_type := fifelse(alteration == "activation", "Activation",
                   fifelse(alteration == "inhibition", "Inhibition", "Bidirectional"))]
e2[, is_core := grepl("CXCL8|PPARG|IL1B|FFAR2", to)]

all_seg <- rbindlist(list(e1[, .(x1, y1, x2, y2, edge_type, is_core)], e2[, .(x1, y1, x2, y2, edge_type, is_core)]))
sig_net <- make_sigmoid_edges(all_seg, direction = "x", n_pts = 70)

node_pos[, node_fill := fifelse(grepl("CXCL8", node), "#FFE4E6",
                        fifelse(grepl("PPARG", node), "#DBEAFE",
                        fifelse(grepl("Excluded|Excl\\.", node), "#F1F5F9",
                        fifelse(col == 1, "#F0F9FF",
                        fifelse(col == 2, "#FEF9C3", "#ECFDF5")))))]
node_pos[, node_bcol := fifelse(grepl("CXCL8", node), pal_npg[["red"]],
                        fifelse(grepl("PPARG", node), pal_npg[["blue"]],
                        fifelse(grepl("Excluded|Excl\\.", node), "#CBD5E1", "#CBD5E1")))]
node_pos[, node_tcol := fifelse(grepl("CXCL8", node), "#9F1239",
                        fifelse(grepl("PPARG", node), "#1E3A8A",
                        fifelse(grepl("Excluded|Excl\\.", node), "#94A3B8", "#0F172A")))]

p2 <- ggplot() +
  geom_path(data = sig_net[is_core == FALSE], aes(x = x, y = y, group = edge_id, color = edge_type, linetype = edge_type),
            linewidth = 0.60, alpha = 0.40) +
  geom_path(data = sig_net[is_core == TRUE],  aes(x = x, y = y, group = edge_id, color = edge_type, linetype = edge_type),
            linewidth = 1.05, alpha = 0.88)

for (i in seq_len(nrow(node_pos))) {
  r <- node_pos[i]
  p2 <- add_round_card(p2, r$col - r$hw, r$col + r$hw, r$y - 0.38, r$y + 0.38,
                       fill = r$node_fill, border_col = r$node_bcol, lwd = ifelse(grepl("CXCL8|PPARG", r$node), 1.1, 0.55), r_pt = 9)
}

col_hdr <- data.table(x = c(1.0, 2.0, 3.0), y = 10.6,
                      lbl = c("Microbial / Metabolite Exposure", "Bioactive Gut Mediators", "Host Genetic Targets"))

p2 <- p2 +
  geom_text(data = col_hdr, aes(x = x, y = y, label = lbl), size = 3.6, fontface = "bold", color = "#334155") +
  geom_text(data = node_pos, aes(x = col, y = y, label = node, color = node_tcol),
            size = 3.05, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(
    values = c("Produces" = "#64748B", "Activation" = pal_npg[["green"]], "Inhibition" = pal_npg[["red"]],
               "Bidirectional" = "#7C3AED", "#9F1239" = "#9F1239", "#1E3A8A" = "#1E3A8A", "#94A3B8" = "#94A3B8", "#0F172A" = "#0F172A"),
    breaks = c("Produces", "Activation", "Inhibition", "Bidirectional"), name = "Interaction"
  ) +
  scale_linetype_manual(
    values = c("Produces" = "solid", "Activation" = "solid", "Inhibition" = "dashed", "Bidirectional" = "dotted"),
    breaks = c("Produces", "Activation", "Inhibition", "Bidirectional"), name = "Interaction"
  ) +
  coord_cartesian(xlim = c(0.62, 3.38), ylim = c(0.8, 11.1), expand = FALSE) +
  labs(title = "Subtype-Divergent Microbiota\u2013Metabolite\u2013Host Causal Network") +
  theme_void(base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 12.2, color = "#0F172A", margin = margin(b = 8)),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 9.5, color = "#0F172A"),
    legend.text      = element_text(size = 9.0, color = "#334155"),
    plot.margin      = margin(10, 12, 10, 12)
  )
save_pub(p2, "Fig2_PlanA_DualAxis_Network", w = 15.2, h = 7.8)
rm(p2, all_seg, sig_net, e1, e2, node_pos, col_hdr); gc()

# ==============================================================================
# FIG 3: 顺式 cis-eQTL MR (语义严格区分 CXCL8红 vs PPARG蓝) + pQTL + Coloc + 峰图
# ==============================================================================
cis_pos[, label := sprintf("%s \u2192 %s (%d %s)", gene_symbol, clean_subtype_str(outcome_label), nsnp, ifelse(nsnp == 1, "SNP", "SNPs"))]
cis_pos[, sig_tag := fifelse(gene_symbol == "CXCR2", "Excluded (LD, PP.H3 = 1.0)",
                     fifelse(gene_symbol == "CXCL8", "CXCL8 (ER\u2212 Risk Axis)",
                     fifelse(gene_symbol == "PPARG", "PPARG (ER+ Protective Axis)",
                     fifelse(FDR < 0.05, "IL1B (FDR < 0.05)", "Nominal (P < 0.05)"))))]
cis_pos[, domain_flag := fifelse(gene_symbol == "CXCR2", "LD-Confounded (Excluded)", "Primary")]
cis_pos[, anno_txt := sprintf("OR = %.3f [%.3f\u2013%.3f], P = %.2g", OR, OR_LCI95, OR_UCI95, pval)]
cis_pos <- cis_pos[order(pval)]
cis_pos[, label := factor(label, levels = rev(label))]

p3a <- plot_forest_panel(
  cis_pos, "A. Strict cis-eQTL \u2192 BCAC MR",
  c(0.45, max(cis_pos$OR_UCI95) * 1.56),
  c("CXCL8 (ER\u2212 Risk Axis)" = pal_npg[["red"]],
    "PPARG (ER+ Protective Axis)" = pal_npg[["blue"]],
    "IL1B (FDR < 0.05)" = pal_npg[["green"]],
    "Nominal (P < 0.05)" = pal_npg[["slate"]],
    "Excluded (LD, PP.H3 = 1.0)" = "#CBD5E1")
)

pqtl_dt <- fread(file.path(res_dir, "planA_pqtl_il8_mr.csv"))
pq_sub <- pqtl_dt[method %in% c("Inverse variance weighted", "Wald ratio") &
                  outcome_label %in% c("BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall")]
pq_sub <- pq_sub[
  (iv_scope == "Cis_1Mb_CXCL8" & !(id.exposure == "ebi-a-GCST90011994" & abs(thresh - 1e-5) < 1e-9)) |
  (iv_scope == "All_IVs" & abs(thresh - 1e-5) < 1e-9)
]
stopifnot(nrow(pq_sub) == 21L)

pq_sub[, cohort_lbl := fifelse(id.exposure == "ebi-a-GCST90011994" & iv_scope == "Cis_1Mb_CXCL8", "SCALLOP cis-pQTL (rs7655660)",
                       fifelse(id.exposure == "ebi-a-GCST90011994", "SCALLOP All IVs (31 SNPs)",
                       fifelse(id.exposure == "prot-a-749" & iv_scope == "Cis_1Mb_CXCL8", "INTERVAL cis-pQTL (rs74327761)",
                       fifelse(id.exposure == "prot-a-749", "INTERVAL All IVs (24 SNPs)",
                       fifelse(id.exposure == "ebi-a-GCST004445", "YFS/FINRISK All IVs (13 SNPs)",
                       fifelse(id.exposure == "prot-c-3447_64_2", "KORA All IVs (7 SNPs)",
                               "LBC All IVs (10 SNPs)"))))))]
pq_sub[, tier := fifelse(pval < 0.05, "P < 0.05", "P \u2265 0.05")]
pq_sub[, out_clean := factor(outcome_label, levels = c("BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall"),
                             labels = c("ER\u2212 Breast Cancer", "ER+ Breast Cancer", "Overall Breast Cancer"))]
ord_pq <- c("SCALLOP cis-pQTL (rs7655660)", "INTERVAL cis-pQTL (rs74327761)",
            "SCALLOP All IVs (31 SNPs)", "INTERVAL All IVs (24 SNPs)",
            "YFS/FINRISK All IVs (13 SNPs)", "LBC All IVs (10 SNPs)", "KORA All IVs (7 SNPs)")
pq_sub[, cohort_lbl := factor(cohort_lbl, levels = rev(ord_pq))]

p3b <- ggplot(pq_sub, aes(x = OR, y = cohort_lbl, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.72) +
  geom_point(size = 2.75) +
  geom_text(aes(x = OR_UCI95 + 0.035, label = ifelse(pval < 0.05, sprintf("%.3f (P = %.3f)", OR, pval), "")),
            hjust = 0, size = 2.6, fontface = "bold", color = "#BE123C", show.legend = FALSE) +
  facet_wrap(~ out_clean, ncol = 3) +
  scale_color_manual(values = c("P < 0.05" = pal_npg[["red"]], "P \u2265 0.05" = pal_npg[["null"]])) +
  coord_cartesian(xlim = c(0.84, 2.34)) +
  labs(title = "B. Plasma CXCL8/IL-8 pQTL \u2192 BCAC MR", x = "Odds Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
coloc_dt[, pair := paste0(gene_symbol, " (", clean_subtype_str(outcome_label), ")")]
pq_coloc <- fread(file.path(res_dir, "planA_pqtl_il8_coloc.csv"))
coloc_comb <- rbind(
  coloc_dt[, .(pair, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)],
  data.table(pair = "IL-8 pQTL (SCALLOP ER\u2212)", PP.H0 = pq_coloc$PP.H0.abf[1], PP.H1 = pq_coloc$PP.H1.abf[1],
             PP.H2 = pq_coloc$PP.H2.abf[1], PP.H3 = pq_coloc$PP.H3.abf[1], PP.H4 = pq_coloc$PP.H4.abf[1])
)[order(PP.H4)]
coloc_long <- melt(coloc_comb, id.vars = "pair", measure.vars = c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),
                   variable.name = "Hypothesis", value.name = "PP")
coloc_long[, pair := factor(pair, levels = coloc_comb$pair)]

p3c <- ggplot(coloc_long, aes(x = PP, y = pair, fill = Hypothesis)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("PP.H0"="#E2E8F0","PP.H1"="#CBD5E1","PP.H2"="#64748B","PP.H3"="#F97316","PP.H4"=pal_npg[["green"]])) +
  labs(title = "C. Bayesian Colocalization (\u00b1500 kb)", x = "Posterior Probability", y = NULL, fill = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

reg_snps <- readRDS(file.path(cache_dir, "planA_coloc_plot_snps.rds"))
lead_sum <- fread(file.path(res_dir, "planA_coloc_lead_snps_summary.csv"))

plot_one_locus <- function(pair_key, title_lbl, show_ylab = FALSE) {
  sub <- reg_snps[pair == pair_key]
  ls  <- lead_sum[pair == pair_key]
  sub[, pos_mb := pos / 1e6]
  comb <- rbind(
    data.table(pos_mb = sub$pos_mb, val = sub$log10p_e, track = "eQTLGen", rsid = sub$rsid),
    data.table(pos_mb = sub$pos_mb, val = sub$log10p_o, track = "BCAC GWAS", rsid = sub$rsid)
  )
  comb[, track := factor(track, levels = c("eQTLGen", "BCAC GWAS"))]
  comb[, is_lead := (track == "eQTLGen" & rsid == ls$lead_eqtl_snp) | (track == "BCAC GWAS" & rsid == ls$lead_brca_snp)]

  ggplot(comb, aes(x = pos_mb, y = val)) +
    ggrastr::rasterise(geom_point(data = comb[is_lead == FALSE], aes(color = track), size = 0.65, alpha = 0.45), dpi = 350) +
    geom_point(data = comb[is_lead == TRUE], color = pal_npg[["red"]], shape = 18, size = 3.0) +
    ggrepel::geom_text_repel(data = comb[is_lead == TRUE], aes(label = rsid), size = 2.6, fontface = "bold", color = "#991B1B",
                             max.overlaps = 10, box.padding = 0.65, point.padding = 0.35, force = 6, force_pull = 0.15,
                             min.segment.length = 0, seed = 123) +
    facet_grid(track ~ ., scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.34))) +
    scale_color_manual(values = c("eQTLGen" = pal_npg[["blue"]], "BCAC GWAS" = pal_npg[["green"]]), guide = "none") +
    labs(title = title_lbl,
         x = sprintf("Chr%s Position (Mb)", sub$chr[1]),
         y = if (show_ylab) expression(bold(-log[10](italic(P)))) else NULL) +
    theme_panel_grid(9.0)
}

p3d <- patchwork::wrap_elements(full = wrap_plots(
  plot_one_locus("PPARG_ER_pos",  "D1. PPARG Locus (ER+)",  show_ylab = TRUE),
  plot_one_locus("IL1B_ER_neg",   "D2. IL1B Locus (ER\u2212)",  show_ylab = FALSE),
  plot_one_locus("CXCL8_ER_neg",  "D3. CXCL8 Locus (ER\u2212)", show_ylab = FALSE),
  plot_one_locus("CXCR2_Overall", "D4. CXCR2 Locus (Excl.)", show_ylab = FALSE),
  ncol = 4
))

fig3 <- (p3a | p3c) / p3b / p3d + plot_layout(heights = c(1.05, 0.90, 1.12))
save_pub(fig3, "Fig3_CisMR_pQTL_Coloc_and_Peaks", w = 17.5, h = 15.0)
rm(fig3, p3a, p3b, p3c, p3d, reg_snps); gc()

# ==============================================================================
# FIG 4: 临床表达分化与跨队列预后 (开放L型主题 + 斜体P + 风险表大间距)
# ==============================================================================
tc  <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb  <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20 <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42 <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]

df_tn <- rbind(
  data.table(Cohort = "GSE42568",
             G_simple = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal", "Tumor"),
             Expr = as.numeric(g42$expr_target["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Normal", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Tumor",  Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)
p4a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple)) +
  geom_violin(alpha = 0.35, trim = FALSE, color = "#475569", linewidth = 0.4) +
  geom_boxplot(width = 0.28, outlier.size = 0.45, alpha = 0.88) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#38BDF8", "Tumor" = pal_npg[["red"]])) +
  labs(title = "A. PPARG Expression (Tumor vs. Normal)", x = NULL, y = expression(bold(Normalized~log[2]~Expression))) +
  theme_journal(9.4) + theme(legend.position = "none")

pam_tc <- rbind(
  data.table(Gene = "CXCL8", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ])),
  data.table(Gene = "PPARG", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)[Subtype %in% c("Normal","LumA","LumB","Her2","Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal","LumA","LumB","Her2","Basal"),
                           labels = c("Normal","LumA","LumB","HER2+","Basal"))]
p4b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_boxplot(width = 0.55, outlier.size = 0.45, alpha = 0.88) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal"="#94A3B8","LumA"=pal_npg[["blue"]],"LumB"="#6366F1","HER2+"=pal_npg[["amber"]],"Basal"=pal_npg[["red"]])) +
  labs(title = "B. PAM50 Subtype Expression (TCGA-BRCA)", x = NULL, y = expression(bold(log[2](RSEM + 1)))) +
  theme_journal(9.4) + theme(legend.position = "none")

idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp],
                      expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4c <- plot_km_with_risk(
  df_tc_p, "time", "status", "Group",
  "C. TCGA-BRCA ER+ OS (PPARG)", "Months", "Overall Survival", "PPARG",
  c("Low", "High"), c(pal_npg[["red"]], pal_npg[["blue"]]),
  pval_coord = c(8, 0.14), break_by = 50, max_x = 250
)

idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp],
                      expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4d <- plot_km_with_risk(
  df_mb_p, "time", "status", "Group",
  "D. METABRIC ER+ OS (PPARG)", "Months", "Overall Survival", "PPARG",
  c("Low", "High"), c(pal_npg[["red"]], pal_npg[["blue"]]),
  pval_coord = c(10, 0.14), break_by = 100, max_x = 300
)

pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4e <- plot_km_with_risk(
  df_g20_c, "time", "status", "Group",
  "E. GSE20685 OS (CXCL8)", "Years", "Overall Survival", "CXCL8",
  c("Low", "High"), c(pal_npg[["blue"]], pal_npg[["red"]]),
  pval_coord = c(0.5, 0.16), break_by = 3, max_x = 15
)

idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern],
                      expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_4f <- plot_km_with_risk(
  df_mb_c, "time", "status", "Group",
  "F. METABRIC ER\u2212 OS (CXCL8)", "Months", "Overall Survival", "CXCL8",
  c("Low", "High"), c(pal_npg[["blue"]], pal_npg[["red"]]),
  pval_coord = c(10, 0.14), break_by = 100, max_x = 300
)
rm(tc, mb, g20, g42, cl_tc); gc()

meta_hr <- fread(file.path(res_dir, "planA_multicohort_meta_survival.csv"))
het_dt  <- fread(file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))
hr_sub <- rbind(
  meta_hr[analysis == "PPARG_5Cohort_Overall_OS", .(label = sprintf("PPARG Overall OS (5 Cohorts, I\u00b2 = %.0f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "ER+ Protective Meta")],
  meta_hr[analysis == "PPARG_2Cohort_ERpos_OS",  .(label = sprintf("PPARG ER+ OS (2 Cohorts, I\u00b2 = %.0f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "ER+ Protective Meta")],
  het_dt[analysis == "All_5_Cohorts_Overall_OS", .(label = sprintf("CXCL8 Overall OS (All 5 Cohorts, I\u00b2 = %.1f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Heterogeneous")],
  het_dt[analysis == "LOO_Excluding_TCGA_BRCA",  .(label = sprintf("CXCL8 Microarray OS (4 Cohorts, I\u00b2 = %.0f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "ER\u2212 Risk Meta")],
  het_dt[analysis == "GEO_3Cohorts_Affy_OS",     .(label = sprintf("CXCL8 Affymetrix GEO OS (3 Cohorts, I\u00b2 = %.0f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "Zero-Het ER\u2212 Meta")],
  meta_hr[analysis == "CXCL8_5Cohort_Progression_RFS_MFS_PFS", .(label = sprintf("CXCL8 Progression RFS/MFS/PFS (5 Cohorts, I\u00b2 = %.1f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "ER\u2212 Risk Meta")],
  meta_hr[analysis == "CXCL8_4Ext_RFS_MFS",      .(label = sprintf("CXCL8 External RFS/MFS (4 Cohorts, I\u00b2 = %.1f%%)", I2_pct), HR = RE_HR, HR_LCI95 = RE_LCI, HR_UCI95 = RE_UCI, pval = RE_pval, tier = "ER\u2212 Risk Meta")]
)
hr_sub[, ord := factor(label, levels = rev(label))]
p4g <- ggplot(hr_sub, aes(x = HR, y = ord, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.28, linewidth = 0.78) +
  geom_point(size = 3.2, shape = 18) +
  geom_text(aes(x = 1.40, label = sprintf("RE HR = %.3f (%.3f\u2013%.3f), P = %.2g", HR, HR_LCI95, HR_UCI95, pval)),
            hjust = 0, size = 2.75, color = "#0F172A", show.legend = FALSE) +
  scale_color_manual(values = c("Zero-Het ER\u2212 Meta" = pal_npg[["red"]], "ER\u2212 Risk Meta" = "#FB7185",
                                "ER+ Protective Meta" = pal_npg[["blue"]], "Heterogeneous" = pal_npg[["null"]])) +
  coord_cartesian(xlim = c(0.74, 1.88)) +
  labs(title = "G. Cross-Cohort Pooled Survival Meta-Analysis",
       x = "Pooled Random-Effects Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

fig4_r1 <- patchwork::wrap_elements(full = (p4a | p4b))
fig4_all <- fig4_r1 / (km_4c | km_4d) / (km_4e | km_4f) / p4g + plot_layout(heights = c(0.85, 1.25, 1.25, 0.95))
save_pub(fig4_all, "Fig4_Clinical_Expression_and_Prognosis", w = 16.5, h = 19.2)
rm(fig4_all, fig4_r1, p4a, p4b, km_4c, km_4d, km_4e, km_4f, p4g); gc()

# ==============================================================================
# FIG 5: 脂肪去混杂多因素验证 (5A) + TCGA ER+ 列线图 (5B) + 增量 (5C) + 校准 (5D)
# ==============================================================================
ad_dt <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
ad_sub <- ad_dt[cohort %in% c("TCGA_ERpos", "TCGA_ERpos_LowADIPOQ", "TCGA_ERpos_LowPLIN1",
                              "METABRIC_ERpos", "METABRIC_ERpos_LowADIPOQ", "METABRIC_ERpos_LowPLIN1")]
ad_sub[, coh_grp := ifelse(grepl("TCGA", cohort), "TCGA-BRCA ER+ (RNA-seq)", "METABRIC ER+ (Microarray)")]
ad_sub[, lbl := fifelse(model == "M1_Univariate", "Univariate PPARG",
                fifelse(model == "M2_Adj_Adipose(ADIPOQ+PLIN1)", "+ ADIPOQ & PLIN1",
                fifelse(model == "M3_Adj_Adipose+Stromal", "+ Adipose & Stromal Score",
                fifelse(model == "M4_Adj_Age+Stage", "+ Age & Stage",
                fifelse(model == "M3_Adj_Age+Stage", "+ Age & Stage",
                fifelse(model == "M4_Adj_Age+Stage+Grade", "+ Age, Stage & Grade",
                fifelse(grepl("M5", model), "Full Multivariable Model",
                fifelse(model == "M6_LowAdipose_Subgroup(ADIPOQ<med)", "Low-ADIPOQ Stratum",
                        "Low-PLIN1 Stratum"))))))))]
ad_sub[, sig := ifelse(pval < 0.05, "P < 0.05", "Attenuated (P \u2265 0.05)")]
ad_sub[, row_ord := factor(paste0(coh_grp, "___", lbl), levels = rev(paste0(coh_grp, "___", lbl)))]

p5a <- ggplot(ad_sub, aes(x = HR, y = row_ord, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
  geom_errorbar(aes(xmin = HR_LCI95, xmax = HR_UCI95), orientation = "y", width = 0.26, linewidth = 0.74) +
  geom_point(size = 2.75) +
  geom_text(aes(x = 1.23, label = sprintf("%.3f (P = %.3f)", HR, pval)), hjust = 0, size = 2.65, show.legend = FALSE) +
  facet_wrap(~ coh_grp, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_color_manual(values = c("P < 0.05" = pal_npg[["blue"]], "Attenuated (P \u2265 0.05)" = pal_npg[["null"]])) +
  coord_cartesian(xlim = c(0.54, 1.56)) +
  labs(title = "A. Adipose-Deconfounded Multivariable Cox (PPARG ER+ OS)", x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

tc_ad <- readRDS(file.path(bulk_dir, "tcga_planA_adipose_clin.rds"))
d_nom <- as.data.frame(tc_ad[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event) & !is.na(AGE) & !is.na(STAGE) & !is.na(PPARG),
                             .(OS_MONTHS, OS_event, AGE, STAGE = factor(STAGE), PPARG = round(PPARG, 2))])
dd_nom <- datadist(d_nom)
assign("dd_nom", dd_nom, envir = .GlobalEnv)
options(datadist = "dd_nom")
fit_nom  <- cph(Surv(OS_MONTHS, OS_event) ~ AGE + STAGE + PPARG, data = d_nom, x = TRUE, y = TRUE, surv = TRUE)
sv_nom   <- Survival(fit_nom)
nom_plot <- nomogram(fit_nom,
                     fun = list(function(x) sv_nom(36, x), function(x) sv_nom(60, x), function(x) sv_nom(120, x)),
                     funlabel = c("3-Year OS", "5-Year OS", "10-Year OS"), lp = FALSE)
tmp_nom_png <- file.path(pub_dir, "_tmp_nomogram_v4.png")
ragg::agg_png(tmp_nom_png, width = 2400, height = 1350, res = 260, background = "white")
par(mar = c(0.8, 0.5, 0.8, 0.5), family = "sans")
plot(nom_plot, xfrac = 0.22, cex.axis = 0.82, cex.var = 0.90, lmgp = 0.22)
dev.off()
p5b <- img_panel(tmp_nom_png, top_ttl = "B. TCGA-BRCA ER+ Clinical\u2013PPARG Nomogram", trim_bg = TRUE, show_border = FALSE)
unlink(tmp_nom_png)

roc_dt <- fread(file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))[cohort == "TCGA_ERpos" & model != "M1_PPARG_Alone"]
roc_dt[, Model := factor(model,
                         levels = c("M2_Clinical_Base","M3_Clin_Plus_PPARG","M4_Clin_Adipose_Base","M5_Clin_Adipose_Plus_PPARG"),
                         labels = c("Age + Stage","Age + Stage + PPARG","Clinical + Adipose + Stromal","Full Model + PPARG"))]
roc_m <- melt(roc_dt, id.vars = "Model", measure.vars = c("C_index","AUC_3yr","AUC_5yr","AUC_10yr"),
              variable.name = "Metric", value.name = "Score")
roc_m[, Metric := factor(Metric, levels = c("C_index","AUC_3yr","AUC_5yr","AUC_10yr"),
                         labels = c("C-index","3-yr AUC","5-yr AUC","10-yr AUC"))]
p5c <- ggplot(roc_m, aes(x = Metric, y = Score, fill = Model)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", Score)), position = position_dodge(width = 0.78), vjust = -0.35, size = 2.45) +
  scale_fill_manual(values = c("#CBD5E1", pal_npg[["blue"]], "#94A3B8", pal_npg[["green"]])) +
  coord_cartesian(ylim = c(0.60, 0.88)) +
  labs(title = "C. TCGA ER+ Discrimination Gain (C-index & timeROC)", x = NULL, y = "Discrimination Score", fill = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

nom_obj <- readRDS(file.path(bulk_dir, "planA_tcga_erpos_nomogram.rds"))
cal_df <- as.data.table(unclass(nom_obj$cal_5yr))
cal_df[, lci := pmax(0, KM.corrected - 1.96 * std.err)]
cal_df[, uci := pmin(1, KM.corrected + 1.96 * std.err)]
p5d <- ggplot(cal_df, aes(x = mean.predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.7) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.015, color = pal_npg[["blue"]], linewidth = 0.7) +
  geom_line(aes(y = KM, color = "Apparent KM"), linewidth = 0.85, linetype = "dotted") +
  geom_point(aes(y = KM, color = "Apparent KM"), size = 2.5) +
  geom_line(aes(y = KM.corrected, color = "Bootstrap Corrected (B = 200)"), linewidth = 1.05) +
  geom_point(aes(y = KM.corrected, color = "Bootstrap Corrected (B = 200)"), size = 3.1) +
  scale_color_manual(values = c("Apparent KM" = pal_npg[["coral"]], "Bootstrap Corrected (B = 200)" = pal_npg[["blue"]])) +
  coord_fixed(xlim = c(0.42, 1.0), ylim = c(0.42, 1.0)) +
  labs(title = "D. 5-Year OS Bootstrap Calibration", x = "Nomogram-Predicted 5-Yr OS", y = "Observed 5-Yr OS (KM)", color = NULL) +
  theme_journal(9.4) + theme(legend.position = "bottom")

fig5_all <- (p5a | patchwork::wrap_elements(full = p5b)) / patchwork::wrap_elements(full = (p5c | p5d)) + plot_layout(heights = c(1.12, 0.95))
save_pub(fig5_all, "Fig5_AdiposeDeconfounding_and_Nomogram", w = 16.5, h = 12.6)
cat("\nPart 1 (Fig 1 - Fig 5) Complete!\n")
