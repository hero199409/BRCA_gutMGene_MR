# 04_check_pos_taxa_gutmgene.R
library(data.table)

data_dir <- "00_data"
sig_taxa <- c("Adlercreutzia", "Desulfovibrio", "Parabacteroides", "Ruminococcaceae", "Dorea", "Blautia", "Marvinbryantia")
pat <- paste(sig_taxa, collapse = "|")

f1 <- fread(file.path(data_dir, "Gut_Microbe_and_Gene-human.txt"), sep = "\t")
f2 <- fread(file.path(data_dir, "Gut_Microbe_and_Metabolite-human.txt"), sep = "\t")
f3 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
f4 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))
f5 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbe_Metabolite_Host_Gene_PubChem_CID_join.csv"))

cat("=== 1. Direct Microbe -> Host Gene (v1) ===\n")
d1 <- f1[grepl(pat, f1[["Gut Microbiota"]], ignore.case = TRUE)]
print(unique(d1[, c("Gut Microbiota", "Gene", "Alteration")]))

cat("\n=== 2. v1 + v2 Microbe -> Metabolite ===\n")
m_v1 <- f2[grepl(pat, f2[["Gut Microbiota"]], ignore.case = TRUE), .(
  source = "v1",
  microbe = get("Gut Microbiota"),
  metabolite = Metabolite,
  pubchem_cid = as.character(get("Metabolite PubChem CID")),
  hmdb = get("Metabolite HMDB ID")
)]
m_v2 <- f3[grepl(pat, gut_microbiota, ignore.case = TRUE), .(
  source = "v2",
  microbe = gut_microbiota,
  metabolite = metabolite,
  pubchem_cid = as.character(metabolite_pubchem_cid),
  hmdb = metabolite_hmdb
)]
m_all <- unique(rbindlist(list(m_v1, m_v2), use.names = TRUE))
print(m_all)

cat("\n=== 3. v2 Triplet Join (Microbe -> Metabolite -> Host Gene) ===\n")
t_sub <- f5[grepl(pat, gut_microbiota, ignore.case = TRUE)]
cat("Triplet rows:", nrow(t_sub),
    "| unique microbes:", length(unique(t_sub$gut_microbiota)),
    "| unique metabolites:", length(unique(t_sub$microbe_table_metabolite)),
    "| unique host genes:", length(unique(t_sub$host_gene)), "\n")
if (nrow(t_sub) > 0) {
  print(unique(t_sub[, .(gut_microbiota, metabolite = microbe_table_metabolite, host_gene, alteration = metabolite_gene_alteration)]))
}

cat("\n=== 4. Additional Metabolite -> Gene matches from v2 browse (by name/PubChem) ===\n")
m_names <- tolower(trimws(m_all$metabolite))
m_cids <- m_all$pubchem_cid[!is.na(m_all$pubchem_cid) & nzchar(m_all$pubchem_cid)]
g_add <- f4[tolower(trimws(metabolite)) %in% m_names | as.character(metabolite_pubchem_cid) %in% m_cids]
cat("Extra Metabolite->Gene rows:", nrow(g_add), "| unique genes:", length(unique(g_add$gene)), "\n")
if (nrow(g_add) > 0) {
  print(unique(g_add[, .(metabolite, gene, alteration)]))
}
