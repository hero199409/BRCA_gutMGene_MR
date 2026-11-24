# Frontiers Fig 1–2. 仅使用本地结果表和 gutMGene 配对记录。
source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/25_0_cns_plot_helpers.R"))

read_result <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))
nominal_mib <- read_result("mibiogen_bcac_mr_significant_p005")[gwas_name != "Marvinbryantia"]
nominal_met <- read_result("metabolite_bcac_mr_significant_p005")[
  id.exposure %chin% c("met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-446")]
cis_audit <- read_result("module3_all41_genes_iv_cis_audit")
cis_primary <- read_result("module3_eqtlgen_cis_only_mr_primary_full_denominator")
cis_nominal <- cis_primary[is.finite(pval) & pval < 0.05]
with_cis <- uniqueN(cis_audit[is_cis_1mb == TRUE, id.exposure])
with_mr <- uniqueN(cis_primary$id.exposure)
stopifnot(with_cis == 40L, with_mr == 36L, nrow(cis_primary) == 108L,
          !anyDuplicated(cis_primary[, .(id.exposure, outcome_label)]),
          all(cis_primary[, .N, by = outcome_label]$N == 36L))

# A. 工作流中的三条分析路径。数字由重跑分母和坐标审计表分别给出。
workflow <- data.table(
  x = 1:3,
  heading = c("Gut microbiota", "Plasma metabolites", "Host expression"),
  detail = c(
    "MiBioGen screen\nReverse MR check\nCross-cohort comparison",
    "Plasma metabolite screen\nReverse MR check\nPrioritized associations",
    sprintf("%d genes with cis IVs\n%d genes in BCAC cis-MR\n%d gene–outcome tests", with_cis, with_mr, nrow(cis_primary))
  )
)
p1a <- ggplot(workflow, aes(x = x)) +
  geom_tile(aes(y = 0.5), width = 0.95, height = 0.94, fill = "#F1F5F9") +
  geom_text(aes(y = 0.79, label = heading), size = 9 / ggplot2::.pt, fontface = "bold", color = "#0F172A") +
  geom_text(aes(y = 0.41, label = detail), size = 8 / ggplot2::.pt, lineheight = 1.13, color = "#334155") +
  coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0, 1), expand = FALSE) +
  labs(title = "A. Study design and cis-eQTL quality control", x = NULL, y = NULL,
       caption = "cis IVs: within ±1 Mb of the GRCh37 gene body. BCAC outcomes: overall, ER+ and ER−.") +
  theme_frontiers_cns() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

prepare_upstream <- function(d, label_col) {
  d <- copy(d)[order(outcome_label, pval)]
  d[, label := sprintf("%s (%s, %d SNPs)", get(label_col), clean_subtype_str(outcome_label), nsnp)]
  d[, tier := fcase(FDR < 0.05, "FDR < 0.05", FDR < 0.10, "0.05 ≤ FDR < 0.10", default = "Nominal P < 0.05")]
  d[, instrument := fifelse(nsnp == 2L, "2 SNPs", "Other instrument counts")]
  d[, annotation := sprintf("%.2f [%.2f, %.2f]   P = %.2g", OR, OR_LCI95, OR_UCI95, pval)]
  d
}
nominal_met[, short_name := fcase(id.exposure == "met-a-385", "3-Phenylpropionate",
                                  id.exposure == "met-c-893", "Linoleic acid (18:2)", default = trait_name)]
up_colors <- c("FDR < 0.05" = "#087F6D", "0.05 ≤ FDR < 0.10" = "#A54A09", "Nominal P < 0.05" = "#3C5488")
up_shapes <- c("Other instrument counts" = 16, "2 SNPs" = 17)
p1b <- forest_frontiers_cns(prepare_upstream(nominal_mib, "gwas_name"), "B. Gut microbiota and BCAC",
                           "MiBioGen results retained after the reverse MR check", up_colors, up_shapes)
p1c <- forest_frontiers_cns(prepare_upstream(nominal_met, "short_name"), "C. Plasma metabolites and BCAC",
                           "Prioritized associations with nominal P < 0.05", up_colors, up_shapes[1]) + guides(shape = "none")
save_frontiers_cns(p1a / p1b / p1c + plot_layout(heights = c(0.92, 2.42, 1.88)),
                   "Fig1_Design_QC_and_Upstream_MR", height_mm = 191)

# Fig 2. 重跑 MR 和数据库直接配对各自展示。数据库不跨表拼接成路径。
meta <- read_result("module3_microbiome_cross_cohort_meta_strict")[analysis_scope == "same_rank"]
meta_keep <- data.table(core_taxon = c("Dorea", "Adlercreutzia"),
                        outcome_label = c("BRCA_ER_neg", "BRCA_ER_pos"))
meta <- meta[meta_keep, on = .(core_taxon, outcome_label), nomatch = 0]
stopifnot(nrow(meta) == nrow(meta_keep))
meta[, `:=`(OR = RE_OR, OR_LCI95 = RE_LCI, OR_UCI95 = RE_UCI,
            label = sprintf("%s (%s)", core_taxon, clean_subtype_str(outcome_label)),
            tier = fifelse(RE_FDR < 0.05, "FDR < 0.05", "FDR ≥ 0.05"), instrument = "Meta-analysis",
            annotation = sprintf("%.3f [%.3f, %.3f]   P = %.3g", RE_OR, RE_LCI, RE_UCI, RE_pval))]
p2a <- forest_frontiers_cns(meta, "A. Selected cross-cohort microbiota associations",
                          "Same-rank random-effects MR meta-analysis", c(up_colors[1], "FDR ≥ 0.05" = "#3C5488"),
                          c("Meta-analysis" = 16)) + guides(shape = "none")

host <- cis_nominal[(gene_symbol %chin% c("CXCL8", "IL1B", "FFAR2") & outcome_label == "BRCA_ER_neg") |
                   (gene_symbol == "PPARG" & outcome_label == "BRCA_ER_pos")][order(outcome_label, pval)]
stopifnot(nrow(host) == 4L)
host[, `:=`(label = sprintf("%s (%s, %d SNP%s)", gene_symbol, clean_subtype_str(outcome_label), nsnp,
                            ifelse(nsnp == 1L, "", "s")),
            tier = fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05"),
            instrument = fifelse(nsnp == 1L, "1 SNP", "2 SNPs"),
            annotation = sprintf("%.3f [%.3f, %.3f]   P = %.3g", OR, OR_LCI95, OR_UCI95, pval))]
p2b <- forest_frontiers_cns(host, "B. Selected host expression associations",
                          "Strict cis-eQTL MR using the full BCAC denominator", up_colors,
                          c("1 SNP" = 16, "2 SNPs" = 17))

gut_dir <- dirname(cache_dir)
microbe <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
gene <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))
# 保留两个源表中各自确切匹配的名称及原始 PMID，不连接 MR 的暴露与细胞靶点。
make_evidence <- function(d, row_key, col_key, row_names, col_names, title) {
  ev <- d[get(row_key) %chin% row_names & get(col_key) %chin% col_names & !is.na(pmid),
          .(row = get(row_key), col = get(col_key), pmid)]
  ev <- ev[, .(evidence = paste(sort(unique(pmid)), collapse = "\n")), by = .(row, col)]
  ev <- ev[CJ(row = row_names, col = col_names, sorted = FALSE), on = .(row, col)]
  ev[, `:=`(row = factor(row, levels = rev(row_names)), col = factor(col, levels = col_names), found = !is.na(evidence))]
  ggplot(ev, aes(col, row)) +
    geom_tile(aes(fill = found), width = 0.96, height = 0.94, color = "white") +
    geom_text(aes(label = ifelse(found, evidence, "")), size = 8 / ggplot2::.pt,
              color = "#123A3A", lineheight = 0.93) +
    scale_fill_manual(values = c(`TRUE` = "#CAEBE4", `FALSE` = "#F1F5F9"), guide = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_frontiers_cns() + theme(axis.line = element_blank(), axis.ticks = element_blank())
}
met_names <- c("Acetate", "Butyrate", "Succinate")
p2c <- make_evidence(microbe, "gut_microbiota", "metabolite", c("Dorea", "Blautia"), met_names,
                     "C. Taxon–metabolite records")
p2d <- make_evidence(gene, "gene", "metabolite", c("CXCL8", "IL1B", "FFAR2", "PPARG"), met_names,
                     "D. Metabolite–gene records")
fig2 <- p2a / p2b / (p2c | p2d) + plot_layout(heights = c(1, 1.38, 1.16)) +
  plot_annotation(
    title = "Selected MR associations and gutMGene records",
    caption = paste("A–B: point = OR, bar = 95% CI. P values are nominal.",
                    "C–D: filled cells list PMIDs for exact pairs in each gutMGene v2.0 source table.",
                    "Blank cells have no matching row in the selected source table. Pair records do not establish mediation.", sep = "\n"),
    theme = theme(plot.title = element_text(size = 11, face = "bold", margin = margin(b = 6)),
                  plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.16, margin = margin(t = 6)),
                  plot.margin = margin(6, 6, 5, 6)))
save_frontiers_cns(fig2, "Fig2_PlanA_DualAxis_Network", height_mm = 147)
