# 02_fig4_tme_gsea.R
# Figure 4 (was Fig 7A-B3): TCGA TME immune correlations (clean single-line labels, zero overlap)
# and Hallmark GSEA lollipop + running enrichment score curves.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})
suppressPackageStartupMessages({
  library(enrichplot)
})

cat("=== Generating Final Fig 4 (TME Correlations & Hallmark GSEA) ===\n")

tme_map <- c(
  "quanTIseq__B cell"                  = "B cells [quanTIseq]",
  "quanTIseq__Macrophage M1"           = "M1 macrophages [quanTIseq]",
  "MCPcounter__Neutrophil"             = "Neutrophils [MCPcounter]",
  "quanTIseq__NK cell"                 = "NK cells [quanTIseq]",
  "TMB_NONSYNONYMOUS"                  = "Nonsynonymous TMB",
  "MCPcounter__T cell CD8+"            = "CD8+ T cells [MCPcounter]",
  "quanTIseq__Neutrophil"              = "Neutrophils [quanTIseq]",
  "quanTIseq__T cell CD8+"             = "CD8+ T cells [quanTIseq]",
  "MCPcounter__T cell"                 = "T cells [MCPcounter]",
  "MCPcounter__B cell"                 = "B cells [MCPcounter]",
  "MCPcounter__Endothelial cell"       = "Endothelial cells [MCPcounter]",
  "ESTIMATE__stroma score"             = "Stromal score [ESTIMATE]",
  "ESTIMATE__tumor purity"             = "Tumor purity [ESTIMATE]",
  "ESTIMATE__estimate score"           = "ESTIMATE score",
  "MCPcounter__cytotoxicity score"     = "Cytotoxicity [MCPcounter]",
  "quanTIseq__uncharacterized cell"    = "Uncharacterized [quanTIseq]",
  "MCPcounter__Myeloid dendritic cell" = "Myeloid DCs [MCPcounter]"
)

cor5     <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
strata7a <- fread(file.path(res_dir, "module5_tcga_er_sample_counts.csv"))[er_status %chin% c("Negative", "Positive")]
stopifnot(nrow(strata7a) == 2L, !anyDuplicated(strata7a$er_status))
strata7a[, stratum := ifelse(er_status == "Negative", "ER_neg", "ER_pos")]
cor5 <- merge(cor5, strata7a[, .(stratum, n_stratum = n_samples)], by = "stratum")

tme_top <- cor5[((gene_symbol == "CXCL8" & stratum == "ER_neg") | (gene_symbol == "PPARG" & stratum == "ER_pos")) &
                  feature_type %chin% c("TME_Cell_Score", "Checkpoint", "Immune_Checkpoint", "Genomic_Feature") & FDR < 0.05][
                    order(gene_symbol, -abs(rho))][, head(.SD, 10L), by = gene_symbol]
stopifnot(uniqueN(tme_top$gene_symbol) == 2L)

tme_top[, feat_clean := ifelse(feature %chin% names(tme_map), unname(tme_map[feature]), gsub("_", " ", feature))]
tme_top[, facet_lbl := sprintf("%s (n = %d)", ifelse(gene_symbol == "CXCL8", "CXCL8 | ER\u2212", "PPARG | ER+"), n_stratum)]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feat_clean)]
tme_top <- tme_top[order(facet_lbl, rho)]
tme_top[, feat_ord := factor(feat_ord, levels = unique(feat_ord))]
tme_top[, y_local := seq_len(.N), by = facet_lbl]
tme_zebra <- tme_top[y_local %% 2 == 1, .(facet_lbl, y_local)]

p4a <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_rect(data = tme_zebra, aes(ymin = y_local - 0.46, ymax = y_local + 0.46, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = pal_rev$zebra, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = rho, yend = feat_ord), linewidth = 0.78) +
  geom_point(size = 2.2) +
  geom_text(aes(x = ifelse(rho >= 0, rho + 0.04, rho - 0.04),
                hjust = ifelse(rho >= 0, 0, 1),
                label = sprintf("%.2f", rho)),
            size = 2.25, fontface = "bold", color = pal_rev$ink, show.legend = FALSE) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_color_manual(values = c("Positive" = pal_rev$cxcl8, "Negative" = pal_rev$pparg)) +
  coord_cartesian(xlim = c(-0.68, 0.85)) +
  labs(title = "A   Tumor microenvironment correlations (TCGA-BRCA)", x = "Spearman \u03c1", y = NULL) +
  theme_rev(8.5) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7.3),
    panel.spacing = unit(4, "mm")
  )

# Panel B: Hallmark pathway enrichment (lollipop)
gsea_map <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"           = "TNF\u03b1\u2013NF-\u03baB signaling",
  "HALLMARK_MTORC1_SIGNALING"                  = "mTORC1 signaling",
  "HALLMARK_G2M_CHECKPOINT"                    = "G2/M checkpoint",
  "HALLMARK_HYPOXIA"                           = "Hypoxia",
  "HALLMARK_GLYCOLYSIS"                        = "Glycolysis",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Inflammatory response",
  "HALLMARK_E2F_TARGETS"                       = "E2F targets",
  "HALLMARK_ANGIOGENESIS"                      = "Angiogenesis",
  "HALLMARK_ALLOGRAFT_REJECTION"               = "Allograft rejection",
  "HALLMARK_KRAS_SIGNALING_UP"                 = "KRAS signaling up",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "IL-6/JAK/STAT3 signaling",
  "HALLMARK_COAGULATION"                       = "Coagulation",
  "HALLMARK_MYOGENESIS"                        = "Myogenesis"
)

gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene == "CXCL8" & stratum == "ER_neg") | (target_gene == "PPARG" & stratum == "ER_pos")) & p.adjust < 0.05][
  order(target_gene, -abs(NES))][, .SD[1:min(8, .N)], by = target_gene]
gs_top[, clean_id := ifelse(ID %in% names(gsea_map), unname(gsea_map[ID]), gsub("_", " ", sub("^HALLMARK_", "", ID)))]
gs_top[, facet_lbl := ifelse(target_gene == "CXCL8", "CXCL8 | ER\u2212", "PPARG | ER+")]
gs_top[, axis_col := ifelse(target_gene == "CXCL8", "CXCL8", "PPARG")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]
gs_top[, id_ord := factor(id_ord, levels = id_ord)]
gs_top[, y_local := seq_len(.N), by = facet_lbl]
gs_zebra <- gs_top[y_local %% 2 == 1, .(facet_lbl, y_local)]

p4b <- ggplot(gs_top, aes(x = NES, y = id_ord, color = axis_col)) +
  geom_rect(data = gs_zebra, aes(ymin = y_local - 0.45, ymax = y_local + 0.45, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = pal_rev$zebra, color = NA) +
  geom_vline(xintercept = 0, color = pal_rev$muted, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, yend = id_ord), linewidth = 1.0, show.legend = FALSE) +
  geom_point(aes(size = -log10(p.adjust))) +
  geom_text(aes(x = NES + 0.16, label = sprintf("%.2f", NES)), hjust = 0, size = 2.35, fontface = "bold",
            color = pal_rev$ink, show.legend = FALSE) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("CXCL8" = pal_rev$cxcl8, "PPARG" = pal_rev$pparg), guide = "none") +
  scale_size(range = c(2.2, 4.6), name = expression(bold(-log[10](FDR)))) +
  coord_cartesian(xlim = c(0, max(gs_top$NES) + 0.65)) +
  labs(title = "B   MSigDB Hallmark enrichment", x = "Normalized enrichment score (NES)", y = NULL) +
  theme_rev(8.5) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(-3, 0, 0, 0),
    axis.text.y = element_text(size = 7.4),
    panel.spacing = unit(4, "mm")
  )

# Panel C: Running enrichment score curves
gsea_objs <- readRDS(file.path(bulk_dir, "planA_cns_gsea_objects.rds"))
curve_specs <- list(
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", lbl = "CXCL8: EMT",                 col = pal_rev$cxcl8),
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",           lbl = "CXCL8: TNF\u03b1\u2013NF-\u03baB", col = "#D97706"),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ALLOGRAFT_REJECTION",               lbl = "PPARG: Allograft rejection", col = pal_rev$pparg),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ANGIOGENESIS",                      lbl = "PPARG: Angiogenesis",        col = pal_rev$teal)
)
curve_dt_list <- vector("list", length(curve_specs))
tick_dt_list  <- vector("list", length(curve_specs))
for (k in seq_along(curve_specs)) {
  sp <- curve_specs[[k]]
  gi <- as.data.table(enrichplot:::gsInfo(sp$obj, geneSetID = sp$id))
  idx_keep <- unique(c(seq(1, nrow(gi), by = max(1L, floor(nrow(gi) / 420))), which.max(abs(gi$runningScore)), nrow(gi)))
  curve_dt_list[[k]] <- gi[idx_keep, .(x, runningScore, pathway = sp$lbl)]
  tick_dt_list[[k]]  <- gi[position == 1, .(x, pathway = sp$lbl, y0 = -0.05 - (k - 1) * 0.045, y1 = -0.085 - (k - 1) * 0.045)]
}
curve_dt <- rbindlist(curve_dt_list)
tick_dt  <- rbindlist(tick_dt_list)
pw_pal   <- setNames(sapply(curve_specs, `[[`, "col"), sapply(curve_specs, `[[`, "lbl"))

p4c <- ggplot(curve_dt, aes(x = x, y = runningScore, color = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
  geom_line(linewidth = 0.90) +
  geom_segment(data = tick_dt, aes(x = x, xend = x, y = y0, yend = y1, color = pathway), linewidth = 0.25, alpha = 0.75) +
  scale_color_manual(values = pw_pal, name = NULL) +
  labs(title = "C   Representative GSEA running enrichment curves", x = "Gene rank", y = "Running ES") +
  theme_rev(8.5) +
  theme(legend.position = "right", legend.text = element_text(size = 7.6))

fig4_all <- wrap_plots(
  wrap_elements(full = p4a),
  wrap_elements(full = p4b),
  wrap_elements(full = p4c),
  ncol = 1,
  heights = c(72, 62, 46)
)
save_rev_figure(fig4_all, "Fig4_TME_Correlations_and_Hallmark_GSEA", width_mm = 180, height_mm = 182, dpi_png = 350)
cat("Fig 4 complete.\n")
