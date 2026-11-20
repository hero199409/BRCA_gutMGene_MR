library(data.table)

f1 <- fread("00_data/Gut_Microbe_and_Gene-human.txt", sep = "\t")
f2 <- fread("00_data/Gut_Microbe_and_Metabolite-human.txt", sep = "\t")
f3 <- fread("00_data/gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv")

cat("1. Gut_Microbe_and_Gene-human.txt:\n")
cat("   rows =", nrow(f1), "| unique_microbe =", length(unique(f1[["Gut Microbiota"]])), "| unique_gene =", length(unique(f1[["Gene"]])), "\n")
cat("   ranks:", paste(names(table(f1$Classification)), table(f1$Classification), sep="=", collapse=", "), "\n")

cat("2. Gut_Microbe_and_Metabolite-human.txt:\n")
cat("   rows =", nrow(f2), "| unique_microbe =", length(unique(f2[["Gut Microbiota"]])), "| unique_metabolite =", length(unique(f2[["Metabolite"]])), "\n")
cat("   ranks:", paste(names(table(f2$Classification)), table(f2$Classification), sep="=", collapse=", "), "\n")

cat("3. gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv:\n")
cat("   rows =", nrow(f3), "\n")
cat("   nonempty gut_microbiota =", sum(!is.na(f3$gut_microbiota) & nzchar(f3$gut_microbiota)), "\n")
cat("   nonempty metabolite =", sum(!is.na(f3$metabolite) & nzchar(f3$metabolite)), "\n")
cat("   nonempty gene =", sum(!is.na(f3$gene) & nzchar(f3$gene)), "\n")
cat("   unique metabolite =", length(unique(f3$metabolite[!is.na(f3$metabolite) & nzchar(f3$metabolite)])), "\n")
cat("   unique gene =", length(unique(f3$gene[!is.na(f3$gene) & nzchar(f3$gene)])), "\n")
cat("   alterations:", paste(names(table(f3$alteration)), table(f3$alteration), sep="=", collapse=", "), "\n")
