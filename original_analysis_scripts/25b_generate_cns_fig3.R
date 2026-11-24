# Frontiers Fig 3. 统计面板使用重跑结果，区域图使用相同的 GRCh37 原始缓存。
source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/25_0_cns_plot_helpers.R"))
read_result <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))
subtype_order <- c("Overall", "ER+", "ER−")

# A. 重跑分母内所有 nominal 结果。CXCR2 保留为筛查结果，C 面板显示区域后验。
cis <- read_result("module3_eqtlgen_cis_only_mr_primary_full_denominator")
stopifnot(nrow(cis) == 108L, uniqueN(cis$id.exposure) == 36L,
          !anyDuplicated(cis[, .(id.exposure, outcome_label)]))
cis <- cis[is.finite(pval) & pval < 0.05][order(gene_symbol, outcome_label)]
cis[, `:=`(label = sprintf("%s (%s, %d SNP%s)", gene_symbol, clean_subtype_str(outcome_label), nsnp,
                            ifelse(nsnp == 1L, "", "s")),
            tier = fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05"),
            instrument = fifelse(nsnp == 1L, "1 SNP", "≥2 SNPs"),
            annotation = sprintf("%.3f [%.3f, %.3f]   P = %.3g", OR, OR_LCI95, OR_UCI95, pval))]
p3a <- forest_frontiers_cns(cis, "A. Strict cis-eQTL Mendelian randomization",
                          "Nominal results from 108 BCAC gene–outcome tests", 
                          c("FDR < 0.05" = "#087F6D", "Nominal P < 0.05" = "#3C5488"),
                          c("1 SNP" = 16, "≥2 SNPs" = 17))

# B. 预先选定的阈值对应重跑 pQTL 表。实际纳入 SNP 数随结局从 nsnp 读取。
pq <- read_result("planA_pqtl_il8_mr")
pq_audit <- read_result("planA_pqtl_il8_iv_audit")
selection <- data.table(
  id.exposure = c("ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST004445", "ebi-a-GCST90000480", "prot-c-3447_64_2"),
  study = c("SCALLOP", "INTERVAL", "SCALLOP", "INTERVAL", "YFS/FINRISK", "LBC", "KORA"),
  iv_scope = c(rep("Cis_1Mb_CXCL8", 2), rep("All_IVs", 5)),
  thresh = c(5e-8, rep(1e-5, 6)))
pq <- pq[method %chin% c("Wald ratio", "Inverse variance weighted") & outcome_label %chin% unique(cis$outcome_label)]
pq <- pq[selection, on = .(id.exposure, iv_scope, thresh), nomatch = 0]
stopifnot(nrow(pq) == 3L * nrow(selection), !anyDuplicated(pq[, .(id.exposure, iv_scope, outcome_label)]))
audit_single <- unique(pq_audit[is_cis_GRCh37 == TRUE, .(id.exposure, thresh = thresh_used, SNP)])
audit_single <- audit_single[selection[iv_scope == "Cis_1Mb_CXCL8"], on = .(id.exposure, thresh), nomatch = 0]
stopifnot(nrow(audit_single) == 2L, !anyDuplicated(audit_single$id.exposure), all(pq[iv_scope == "Cis_1Mb_CXCL8", nsnp] == 1L))
pq[, row_key := paste(study, iv_scope)]
sel_labels <- selection[, paste(study, iv_scope)]
pq[, row_key := factor(row_key, levels = rev(sel_labels), labels = rev(c(
  paste0("SCALLOP cis (", audit_single[id.exposure == "ebi-a-GCST90011994", SNP], ")"),
  paste0("INTERVAL cis (", audit_single[id.exposure == "prot-a-749", SNP], ")"),
  "SCALLOP all IVs", "INTERVAL all IVs", "YFS/FINRISK all IVs", "LBC all IVs", "KORA all IVs")))]
pq[, `:=`(outcome = factor(clean_subtype_str(outcome_label), levels = subtype_order),
            tier = fifelse(pval < 0.05, "P < 0.05", "P ≥ 0.05"))]
p3b <- ggplot(pq, aes(OR, row_key, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.2, linewidth = 0.55) +
  geom_point(size = 2.2) +
  geom_text(aes(x = exp(max(log(OR_UCI95)) + 0.11), label = nsnp), hjust = 0,
            color = "#334155", size = 8 / ggplot2::.pt, show.legend = FALSE) +
  facet_wrap(vars(outcome), nrow = 1) +
  scale_x_log10(breaks = c(0.8, 1, 1.5, 2), labels = label_number()) +
  scale_color_manual(values = c("P < 0.05" = "#B44432", "P ≥ 0.05" = "#475569")) +
  labs(title = "B. Plasma IL-8 pQTL Mendelian randomization",
       subtitle = "Points and bars show OR and 95% CI. Numbers at right give SNP counts.",
       x = "Odds ratio (log scale)", y = NULL, color = NULL) + theme_frontiers_cns()

# C. 只读重跑的 eQTL 后验及 pQTL 可用性报告。
coloc <- read_result("module3_cis_positive_genes_coloc")[p12_prior == 1e-5]
stopifnot(nrow(coloc) > 0L, all(is.finite(coloc$PP.H3)), all(is.finite(coloc$PP.H4)),
          !anyDuplicated(coloc[, .(gene_symbol, outcome_label)]))
sc_status <- read_result("planA_pqtl_il8_availability")[analysis == "SCALLOP coloc", status]
stopifnot(identical(sc_status, "GRCh37 summary statistics not cached"))
coloc[, pair := sprintf("%s (%s)", gene_symbol, clean_subtype_str(outcome_label))]
coloc <- coloc[order(PP.H4 - PP.H3)]
coloc[, pair := factor(pair, levels = pair)]
mirror <- rbind(coloc[, .(pair, probability = -PP.H3, hypothesis = "PP.H3: distinct signals")],
                coloc[, .(pair, probability =  PP.H4, hypothesis = "PP.H4: shared signal")])
p3c <- ggplot(mirror, aes(probability, pair, color = hypothesis)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "#475569") +
  geom_segment(aes(x = 0, xend = probability, yend = pair), linewidth = 0.55) +
  geom_point(size = 1.9) +
  geom_text(data = mirror[abs(probability) >= 0.05],
            aes(x = probability + sign(probability) * 0.035, hjust = ifelse(probability < 0, 1, 0),
                label = ifelse(abs(probability) > 0.99, ">0.99", sprintf("%.2f", abs(probability)))),
            size = 8 / ggplot2::.pt, show.legend = FALSE) +
  scale_color_manual(values = c("PP.H3: distinct signals" = "#A54A09", "PP.H4: shared signal" = "#087F6D")) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), labels = function(x) format(abs(x), trim = TRUE)) +
  coord_cartesian(xlim = c(-1.15, 1.15)) +
  labs(title = "C. cis-eQTL and BCAC regional colocalization",
       subtitle = "GRCh37 windows of ±500 kb. Shared-signal prior p12 = 10⁻⁵.",
       caption = "SCALLOP IL-8 colocalization is not reported because GRCh37 regional summary statistics are not cached.",
       x = "Posterior probability (mirrored)", y = NULL, color = NULL) + theme_frontiers_cns()

# D. 从 C 面板对应的区域缓存重建共有 SNP 集合，核验与重跑表的 SNP 数完全一致。
# 离线注释缓存仅用于基因/外显子位置。此处使用两种信号的 lead 标线，不读取历史后验或峰距。
annotation <- readRDS(file.path(cache_dir, "planA_cns_loci_anno_cache.rds"))
first_letter <- c("D1", "D2", "D3", "D4")
region_jobs <- data.table(gene = c("PPARG", "IL1B", "CXCL8", "CXCR2"),
                          outcome = c("BRCA_ER_pos", "BRCA_ER_neg", "BRCA_ER_neg", "BRCA_Overall"),
                          id = c("ieu-a-1127", "ieu-a-1128", "ieu-a-1128", "ieu-a-1126"),
                          title = first_letter)
track_colors <- c(eQTLGen = "#3C5488", BCAC = "#B44432")
plot_region <- function(job) {
  row <- coloc[gene_symbol == job$gene & outcome_label == job$outcome]
  stopifnot(nrow(row) == 1L)
  coords <- strsplit(row$region, "[:-]")[[1L]]
  xlim <- as.numeric(coords[2:3]) / 1e6
  input <- readRDS(file.path(cache_dir, sprintf("coloc_%s_%s_%s_%s.rds", row$id.exposure, job$id, coords[1L], coords[2L])))
  clean <- function(x) as.data.table(x)[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
  snps <- merge(clean(input$e_reg), clean(input$o_reg), by = "rsid", suffixes = c(".e", ".o"))[
    is.finite(eaf.e) & eaf.e > 0 & eaf.e < 1]
  stopifnot(nrow(snps) == row$n_shared_snps, all(is.finite(snps$position.e)), all(snps$position.e == snps$position.o))
  tracks <- rbind(snps[, .(snp = rsid, x = position.e / 1e6, p = p.e, study = "eQTLGen")],
                  snps[, .(snp = rsid, x = position.o / 1e6, p = p.o, study = "BCAC")])
  tracks[, study := factor(study, levels = names(track_colors))]
  tracks[, logp := -log10(pmax(p, .Machine$double.xmin))]
  peaks <- tracks[is.finite(p)][order(p, snp), .SD[1], by = study]
  stopifnot(nrow(peaks) == 2L)
  gap <- diff(range(peaks$x)) * 1000
  panel <- ggplot(tracks, aes(x, logp)) +
    geom_vline(data = peaks, aes(xintercept = x, color = study), linetype = "dashed", linewidth = 0.45) +
    geom_point(size = 0.43, color = "#6A7788", alpha = 0.65) +
    geom_point(data = peaks, aes(color = study), size = 1.65) +
    facet_grid(rows = vars(study), scales = "free_y") +
    scale_color_manual(values = track_colors, guide = "none") +
    scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = scales::breaks_extended(3), expand = expansion(mult = c(0.04, 0.16))) +
    labs(title = sprintf("%s. %s (%s)", job$title, job$gene, clean_subtype_str(job$outcome)),
         subtitle = sprintf("Lead SNP distance %.1f kb", gap), x = NULL, y = expression(-log[10](italic(P)))) +
    theme_frontiers_cns() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                  panel.spacing.y = unit(0.55, "mm"), plot.margin = margin(4, 6, 0, 6))
  key <- paste(job$gene, sub("BRCA_", "", job$outcome), sep = "_")
  genes <- as.data.table(annotation$genes)[pair == key & is_target == TRUE]
  exons <- as.data.table(annotation$exons)[pair == key & gene_name == job$gene & start/1e6 <= xlim[2] & end/1e6 >= xlim[1]]
  stopifnot(nrow(genes) == 1L, nrow(exons) > 0L)
  track <- ggplot() +
    geom_segment(data = genes, aes(x = pmax(start/1e6, xlim[1]), xend = pmin(end/1e6, xlim[2]), y = 0.5, yend = 0.5),
                 linewidth = 0.45, color = "#475569") +
    geom_rect(data = exons, aes(xmin = pmax(start/1e6, xlim[1]), xmax = pmin(end/1e6, xlim[2]), ymin = 0.12, ymax = 0.88), fill = "#475569") +
    scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    labs(x = sprintf("Chr%s position (Mb, GRCh37)", coords[1]), y = NULL) +
    theme_frontiers_cns() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                  axis.line.y = element_blank(), plot.margin = margin(0, 6, 2, 6))
  panel / track + plot_layout(heights = c(1, 0.09))
}
regions <- lapply(seq_len(nrow(region_jobs)), function(i) plot_region(region_jobs[i]))
# 在每个区域内对齐两条信号和外显子轨，四个区域采用两列排版。
p3d <- wrap_plots(lapply(regions, function(p) wrap_elements(full = p)), ncol = 2)
fig3 <- p3a / p3b / p3c / p3d + plot_layout(heights = c(2.75, 2.05, 2.82, 4.78)) +
  plot_annotation(caption = paste("D1–D4: dashed lines mark the lead eQTLGen (blue) and BCAC (red) SNPs.",
                                  "The lower track shows target gene exons. Peak distance alone does not establish a shared signal.", sep = "\n"),
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.15, margin = margin(t = 3))))
save_frontiers_cns(fig3, "Fig3_CisMR_pQTL_Coloc_and_Peaks", height_mm = 252)
