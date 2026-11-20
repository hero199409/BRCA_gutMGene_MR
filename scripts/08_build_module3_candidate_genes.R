# 08_build_module3_candidate_genes.R
# 严格汇总所有上游 MR 阳性菌群(6个) + 阳性代谢物(5个) 在 gutMGene 中映射到的完整候选宿主基因分母，并匹配 eQTLGen ID (Ensembl)

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

data_dir <- "00_data"
res_dir  <- "02_results"

f1 <- fread(file.path(data_dir, "Gut_Microbe_and_Gene-human.txt"), sep = "\t")
f2 <- fread(file.path(data_dir, "Gut_Microbe_and_Metabolite-human.txt"), sep = "\t")
f3 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
f4 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))
f5 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbe_Metabolite_Host_Gene_PubChem_CID_join.csv"))
ao <- readRDS(file.path(data_dir, "opengwas_available_outcomes.rds"))

# 6 个有效阳性菌群 (不含反向 MR 阳性的 Marvinbryantia)
valid_taxa <- c("Adlercreutzia", "Desulfovibrio", "Parabacteroides", "Ruminococcaceae", "Dorea", "Blautia")
taxa_pat   <- paste(valid_taxa, collapse = "|")

# 1) 从阳性菌群 -> 代谢物 (v1 + v2)
m_v1 <- f2[grepl(taxa_pat, f2[["Gut Microbiota"]], ignore.case = TRUE), .(
  microbe = get("Gut Microbiota"),
  metabolite = Metabolite,
  pubchem_cid = as.character(get("Metabolite PubChem CID"))
)]
m_v2 <- f3[grepl(taxa_pat, gut_microbiota, ignore.case = TRUE), .(
  microbe = gut_microbiota,
  metabolite = metabolite,
  pubchem_cid = as.character(metabolite_pubchem_cid)
)]
taxa_mets <- unique(rbind(m_v1, m_v2))

# 2) 加入 5 个 MR 阳性血浆代谢物 (含别名: Acetate, Linoleic acid, Indoleacetate/Indole-3-acetic acid 等)
mr_met_names <- c("acetate", "acetic acid", "linoleic acid", "indoleacetate", "indole-3-acetic acid",
                  "3-phenylpropionate", "hydrocinnamate", "3-phenylpropionic acid",
                  "3-(4-hydroxyphenyl)lactate", "4-hydroxyphenyllactic acid")
# Acetate PubChem CID = 175 (or 176 acetic acid), Linoleic acid = 5280450, Indole-3-acetic acid = 802
mr_met_cids <- c("175", "176", "5280450", "802")

all_cand_met_names <- unique(c(tolower(trimws(taxa_mets$metabolite)), mr_met_names))
all_cand_cids <- unique(c(taxa_mets$pubchem_cid[!is.na(taxa_mets$pubchem_cid) & nzchar(taxa_mets$pubchem_cid)], mr_met_cids))

# 3) 提取所有候选调控边:
# 3A: 直接 Microbe -> Gene (v1)
e_direct <- f1[grepl(taxa_pat, f1[["Gut Microbiota"]], ignore.case = TRUE), .(
  source_table = "v1_microbe_gene",
  microbe = get("Gut Microbiota"),
  metabolite = NA_character_,
  gene_symbol = trimws(Gene),
  gene_ncbi_id = as.character(get("Gene ID")),
  alteration = tolower(trimws(Alteration)),
  pmid = NA_character_
)]

# 3B: Triplet 表中命中阳性菌群或阳性代谢物
e_trip <- f5[grepl(taxa_pat, gut_microbiota, ignore.case = TRUE) |
               tolower(trimws(microbe_table_metabolite)) %in% all_cand_met_names |
               tolower(trimws(gene_table_metabolite)) %in% all_cand_met_names |
               as.character(metabolite_pubchem_cid) %in% all_cand_cids, .(
  source_table = "v2_triplet_join",
  microbe = gut_microbiota,
  metabolite = gene_table_metabolite,
  gene_symbol = trimws(host_gene),
  gene_ncbi_id = as.character(host_gene_ncbi_id),
  alteration = tolower(trimws(metabolite_gene_alteration)),
  pmid = as.character(metabolite_gene_pmids)
)]

# 3C: v2 Metabolite -> Gene browse 表中命中上述所有代谢物
e_metgene <- f4[tolower(trimws(metabolite)) %in% all_cand_met_names |
                  as.character(metabolite_pubchem_cid) %in% all_cand_cids, .(
  source_table = "v2_metabolite_gene",
  microbe = NA_character_,
  metabolite = metabolite,
  gene_symbol = trimws(gene),
  gene_ncbi_id = as.character(gene_id),
  alteration = tolower(trimws(alteration)),
  pmid = as.character(pmid)
)]

all_edges <- unique(rbindlist(list(e_direct, e_trip, e_metgene), use.names = TRUE))
all_edges <- all_edges[!is.na(gene_symbol) & nzchar(gene_symbol)]

# 把 v2_metabolite_gene 的上游菌群补注进去（只标注来自 6 个阳性菌群还是 MR 阳性代谢物 Acetate 等）
for (i in seq_len(nrow(all_edges))) {
  if (is.na(all_edges$microbe[i])) {
    mn <- tolower(trimws(all_edges$metabolite[i]))
    hits <- unique(taxa_mets[tolower(trimws(metabolite)) == mn]$microbe)
    if (length(hits) > 0) {
      all_edges$microbe[i] <- paste(unique(sapply(strsplit(hits, " "), `[`, 1)), collapse = ";")
    } else if (mn %in% mr_met_names) {
      all_edges$microbe[i] <- "MR_Positive_Metabolite_Direct"
    }
  }
}

fwrite(all_edges, file.path(res_dir, "module3_candidate_gutmgene_edges.csv"))

# 4) 整理去重的候选宿主基因清单，映射 Ensembl ID 并匹配 eQTLGen
gene_dt <- unique(all_edges[, .(gene_symbol, gene_ncbi_id)])
# 汇总每个基因的上游支撑信息
gene_meta <- all_edges[, .(
  upstream_microbes = paste(unique(na.omit(microbe)), collapse = "; "),
  upstream_metabolites = paste(unique(na.omit(metabolite)), collapse = "; "),
  alterations = paste(unique(na.omit(alteration)), collapse = "; "),
  n_edges = .N
), by = gene_symbol]

gene_dt <- merge(gene_dt, gene_meta, by = "gene_symbol", all.x = TRUE)

# 用 org.Hs.eg.db 把 ENTREZ / SYMBOL 映射到 ENSEMBL
map_entrez <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gene_dt$gene_ncbi_id[!is.na(gene_dt$gene_ncbi_id) & nzchar(gene_dt$gene_ncbi_id)]),
  columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
  keytype = "ENTREZID"
))
map_dt <- as.data.table(map_entrez)

eqtl_ao <- ao[grepl("^eqtl-a-", id)]
eqtl_ao[, ensembl := sub("^eqtl-a-", "", id)]

# 匹配 Ensembl
gene_dt[, eqtl_id := NA_character_]
gene_dt[, ensembl_id := NA_character_]
gene_dt[, eqtl_sample_size := NA_integer_]

for (i in seq_len(nrow(gene_dt))) {
  eid <- gene_dt$gene_ncbi_id[i]
  sym <- gene_dt$gene_symbol[i]
  ens_cands <- unique(na.omit(map_dt[ENTREZID == eid]$ENSEMBL))
  if (length(ens_cands) == 0) {
    m2 <- suppressMessages(tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = sym, columns = "ENSEMBL", keytype = "SYMBOL"), error = function(e) NULL))
    if (!is.null(m2)) ens_cands <- unique(na.omit(m2$ENSEMBL))
  }
  hit <- eqtl_ao[ensembl %in% ens_cands]
  if (nrow(hit) > 0) {
    hit <- hit[order(-sample_size)][1]
    gene_dt$eqtl_id[i] <- hit$id
    gene_dt$ensembl_id[i] <- hit$ensembl
    gene_dt$eqtl_sample_size[i] <- hit$sample_size
  }
}

# 按 gene_symbol 去重保留最优一条
gene_dt <- gene_dt[order(gene_symbol, -eqtl_sample_size)]
gene_dt <- gene_dt[!duplicated(gene_symbol)]

fwrite(gene_dt, file.path(res_dir, "module3_candidate_host_genes_universe.csv"))

cat("=== 模块 3 候选宿主基因分母统计 ===\n")
cat("Total unique candidate host genes from causal microbes/metabolites:", nrow(gene_dt), "\n")
cat("Matched in eQTLGen (eqtl-a-*):", sum(!is.na(gene_dt$eqtl_id)), "\n")
cat("Missing in eQTLGen:", sum(is.na(gene_dt$eqtl_id)), "\n")
if (sum(is.na(gene_dt$eqtl_id)) > 0) {
  cat("Missing genes:", paste(gene_dt[is.na(eqtl_id)]$gene_symbol, collapse = ", "), "\n")
}
