# 02_fig8_cis_mr_coloc.R
# Figure 8 (was Old Fig 3): Blood cis-eQTL MR, plasma IL-8 pQTL MR (fixed right-edge SNP clipping!),
# mirrored regional colocalization, and taller, fully legible GRCh37 locus plots.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Final Fig 8 (cis-MR, pQTL, Coloc & GRCh37 Loci) ===\n")

read_res <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))
subtype_order <- c("Overall", "ER+", "ER\u2212")

# 8A: Strict blood cis-eQTL MR (13 nominal from 108 BCAC tests)
cis <- read_res("module3_eqtlgen_cis_only_mr_primary_full_denominator")
stopifnot(nrow(cis) == 108L, uniqueN(cis$id.exposure) == 36L)
cis <- cis[is.finite(pval) & pval < 0.05][order(gene_symbol, outcome_label)]
cis[, `:=`(
  label = sprintf("%s (%s, %d SNP%s)", gene_symbol, clean_subtype_rev(outcome_label), nsnp, ifelse(nsnp == 1L, "", "s")),
  tier  = fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05"),
  inst  = fifelse(nsnp == 1L, "1 SNP", "\u22652 SNPs"),
  stats = sprintf("%.3f [%.3f, %.3f] | P=%.3g", OR, OR_LCI95, OR_UCI95, pval)
)]
cis[, label := factor(label, levels = rev(label))]
cis[, y_num := as.numeric(label)]
lo_a <- min(log(cis$OR_LCI95), 0); hi_a <- max(log(cis$OR_UCI95), 0); sp_a <- hi_a - lo_a
cis[, ann_x := exp(hi_a + 0.06 * sp_a)]

p8a <- ggplot(cis, aes(x = OR, y = label, color = tier)) +
  geom_rect(data = cis[y_num %% 2L == 1L], aes(ymin = y_num - 0.46, ymax = y_num + 0.46),
            xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_rev$zebra) +
  geom_vline(xintercept = 1, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.24, linewidth = 0.58) +
  geom_point(aes(shape = inst), size = 2.0) +
  geom_text(aes(x = ann_x, label = stats), hjust = 0, size = 2.15, color = pal_rev$ink, show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = pal_rev$teal, "Nominal P < 0.05" = pal_rev$pparg), name = NULL) +
  scale_shape_manual(values = c("1 SNP" = 16, "\u22652 SNPs" = 17), name = NULL) +
  scale_x_log10(breaks = c(0.6, 1.0, 1.5), labels = c("0.6", "1.0", "1.5")) +
  coord_cartesian(xlim = exp(c(lo_a - 0.05 * sp_a, hi_a + 0.98 * sp_a))) +
  labs(title = "A   Blood cis-eQTL MR (108 BCAC tests)", x = "Odds ratio (log scale)", y = NULL) +
  theme_rev(8.1) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        axis.text.y = element_text(size = 7.0))

# 8B: Plasma IL-8 pQTL MR (generous right x-limit so SNP counts NEVER clip!)
pq <- read_res("planA_pqtl_il8_mr")
pq_audit <- read_res("planA_pqtl_il8_iv_audit")
selection <- data.table(
  id.exposure = c("ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST004445", "ebi-a-GCST90000480", "prot-c-3447_64_2"),
  study       = c("SCALLOP", "INTERVAL", "SCALLOP", "INTERVAL", "YFS/FINRISK", "LBC", "KORA"),
  iv_scope    = c(rep("Cis_1Mb_CXCL8", 2), rep("All_IVs", 5)),
  thresh      = c(5e-8, rep(1e-5, 6))
)
pq <- pq[method %chin% c("Wald ratio", "Inverse variance weighted") & outcome_label %chin% c("BRCA_Overall", "BRCA_ER_pos", "BRCA_ER_neg")]
pq <- pq[selection, on = .(id.exposure, iv_scope, thresh), nomatch = 0]
audit_single <- unique(pq_audit[is_cis_GRCh37 == TRUE, .(id.exposure, thresh = thresh_used, SNP)])
audit_single <- audit_single[selection[iv_scope == "Cis_1Mb_CXCL8"], on = .(id.exposure, thresh), nomatch = 0]

pq[, row_key := paste(study, iv_scope)]
sel_labels <- selection[, paste(study, iv_scope)]
pq[, row_key := factor(row_key, levels = rev(sel_labels), labels = rev(c(
  paste0("SCALLOP cis (", audit_single[id.exposure == "ebi-a-GCST90011994", SNP], ")"),
  paste0("INTERVAL cis (", audit_single[id.exposure == "prot-a-749", SNP], ")"),
  "SCALLOP all IVs", "INTERVAL all IVs", "YFS/FINRISK all IVs", "LBC all IVs", "KORA all IVs")))]
pq[, `:=`(outcome = factor(clean_subtype_rev(outcome_label), levels = subtype_order),
          tier    = fifelse(pval < 0.05, "Nominal P < 0.05", "P \u2265 0.05"))]

p8b <- ggplot(pq, aes(OR, row_key, color = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = pal_rev$muted, linewidth = 0.42) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.22, linewidth = 0.58) +
  geom_point(size = 2.0) +
  geom_text(aes(x = 2.18, label = paste0("n=", nsnp)), hjust = 0, size = 2.1, color = pal_rev$subink, show.legend = FALSE) +
  facet_wrap(vars(outcome), nrow = 1) +
  scale_x_log10(breaks = c(0.8, 1.0, 1.5, 2.0), labels = c("0.8", "1.0", "1.5", "2.0")) +
  coord_cartesian(xlim = c(0.75, 2.75)) +
  scale_color_manual(values = c("Nominal P < 0.05" = pal_rev$cxcl8, "P \u2265 0.05" = pal_rev$slate), name = NULL) +
  labs(title = "B   Plasma IL-8 pQTL MR", x = "Odds ratio (log scale)", y = NULL) +
  theme_rev(8.1) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        axis.text.y = element_text(size = 7.0), panel.spacing = unit(3, "mm"))

# 8C: Mirrored colocalization lollipop (p12 = 1e-5)
coloc <- read_res("module3_cis_positive_genes_coloc")[p12_prior == 1e-5]
coloc[, pair := sprintf("%s (%s)", gene_symbol, clean_subtype_rev(outcome_label))]
coloc <- coloc[order(PP.H4 - PP.H3)]
coloc[, pair := factor(pair, levels = pair)]
mirror <- rbind(coloc[, .(pair, probability = -PP.H3, hypothesis = "PP.H3 (distinct)")],
                coloc[, .(pair, probability =  PP.H4, hypothesis = "PP.H4 (shared)")])

p8c <- ggplot(mirror, aes(probability, pair, color = hypothesis)) +
  geom_vline(xintercept = 0, linewidth = 0.42, color = pal_rev$subink) +
  geom_segment(aes(x = 0, xend = probability, yend = pair), linewidth = 0.58) +
  geom_point(size = 1.85) +
  geom_text(data = mirror[abs(probability) >= 0.05],
            aes(x = probability + sign(probability) * 0.04, hjust = ifelse(probability < 0, 1, 0),
                label = ifelse(abs(probability) > 0.99, ">0.99", sprintf("%.2f", abs(probability)))),
            size = 2.1, show.legend = FALSE) +
  scale_color_manual(values = c("PP.H3 (distinct)" = pal_rev$amber, "PP.H4 (shared)" = pal_rev$teal), name = NULL) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), labels = c("1.0", "0.5", "0.0", "0.5", "1.0")) +
  coord_cartesian(xlim = c(-1.18, 1.18)) +
  labs(title = "C   cis-eQTL and BCAC colocalization (\u00b1500 kb, p12 = 10\u207b\u2075)",
       x = "Posterior probability (mirrored)", y = NULL) +
  theme_rev(8.1) +
  theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
        axis.text.y = element_text(size = 7.0))

# 8D: Spacious, tall GRCh37 locus plots (D1-D4; legible y-ticks & strip labels!)
annotation <- readRDS(file.path(cache_dir, "planA_cns_loci_anno_cache.rds"))
region_jobs <- data.table(
  gene    = c("CXCL8", "IL1B", "PPARG", "CXCR2"),
  outcome = c("BRCA_ER_neg", "BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall"),
  id      = c("ieu-a-1128", "ieu-a-1128", "ieu-a-1127", "ieu-a-1126"),
  title   = c("D1", "D2", "D3", "D4")
)
track_colors <- c(eQTLGen = pal_rev$pparg, BCAC = pal_rev$cxcl8)

plot_region_rev <- function(job) {
  row <- coloc[gene_symbol == job$gene & outcome_label == job$outcome]
  coords <- strsplit(row$region, "[:-]")[[1L]]
  xlim <- as.numeric(coords[2:3]) / 1e6
  input <- readRDS(file.path(cache_dir, sprintf("coloc_%s_%s_%s_%s.rds", row$id.exposure, job$id, coords[1L], coords[2L])))
  clean_f <- function(x) as.data.table(x)[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
  snps <- merge(clean_f(input$e_reg), clean_f(input$o_reg), by = "rsid", suffixes = c(".e", ".o"))[
    is.finite(eaf.e) & eaf.e > 0 & eaf.e < 1]
  tracks <- rbind(snps[, .(snp = rsid, x = position.e / 1e6, p = p.e, study = "eQTLGen")],
                  snps[, .(snp = rsid, x = position.o / 1e6, p = p.o, study = "BCAC")])
  tracks[, study := factor(study, levels = names(track_colors))]
  tracks[, logp := -log10(pmax(p, .Machine$double.xmin))]
  peaks <- tracks[is.finite(p)][order(p, snp), .SD[1], by = study]
  gap <- diff(range(peaks$x)) * 1000

  p_top <- ggplot(tracks, aes(x, logp)) +
    geom_vline(data = peaks, aes(xintercept = x, color = study), linetype = "dashed", linewidth = 0.45) +
    ggrastr::rasterise(geom_point(size = 0.45, color = "#64748B", alpha = 0.65), dpi = 300) +
    geom_point(data = peaks, aes(color = study), size = 1.8) +
    facet_grid(rows = vars(study), scales = "free_y") +
    scale_color_manual(values = track_colors, guide = "none") +
    scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = function(lims) { b <- scales::breaks_pretty(n = 2)(c(0, lims[2])); b[b >= 0] },
                       expand = expansion(mult = c(0.06, 0.18))) +
    labs(title = sprintf("%s  %s (%s) | \u0394 = %.1f kb", job$title, job$gene, clean_subtype_rev(job$outcome), gap),
         x = NULL, y = expression(-log[10](italic(P)))) +
    theme_rev(7.6) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 6.5),
      strip.text.y = element_text(size = 6.6, margin = margin(l = 2, r = 2)),
      panel.spacing.y = unit(1.2, "mm"),
      plot.title   = element_text(size = 7.9, margin = margin(b = 2)),
      plot.margin  = margin(2, 4, 0, 4)
    )

  key <- paste(job$gene, sub("BRCA_", "", job$outcome), sep = "_")
  genes <- as.data.table(annotation$genes)[pair == key & is_target == TRUE]
  exons <- as.data.table(annotation$exons)[pair == key & gene_name == job$gene & start / 1e6 <= xlim[2] & end / 1e6 >= xlim[1]]
  p_ex <- ggplot() +
    geom_segment(data = genes, aes(x = pmax(start / 1e6, xlim[1]), xend = pmin(end / 1e6, xlim[2]), y = 0.5, yend = 0.5),
                 linewidth = 0.5, color = pal_rev$subink) +
    geom_rect(data = exons, aes(xmin = pmax(start / 1e6, xlim[1]), xmax = pmin(end / 1e6, xlim[2]), ymin = 0.15, ymax = 0.85),
              fill = pal_rev$subink) +
    scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = sprintf("Chr%s (Mb, GRCh37)", coords[1]), y = NULL) +
    theme_rev(7.6) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_text(size = 6.7), axis.title.x = element_text(size = 7.1),
          plot.margin = margin(0, 4, 2, 4))

  wrap_elements(full = (p_top / p_ex + plot_layout(heights = c(1, 0.14))))
}

regions <- lapply(seq_len(nrow(region_jobs)), function(i) plot_region_rev(region_jobs[i]))
p8d <- wrap_elements(full = wrap_plots(regions, ncol = 2))

fig8_all <- wrap_plots(
  wrap_elements(full = p8a),
  wrap_elements(full = p8b),
  wrap_elements(full = p8c),
  p8d,
  ncol = 1,
  heights = c(56, 44, 52, 92)
)
save_rev_figure(fig8_all, "Fig8_CisMR_pQTL_Coloc_and_Loci", width_mm = 180, height_mm = 244, dpi_png = 350)
cat("Fig 8 complete.\n")
