# 17c_extract_gse176078_targets.py
# Fast, low-RAM (<100 MB) binary chunk extractor for Wu et al. 2021 (GSE176078)
# Reads 2.4 GB MatrixMarket count_matrix_sparse.mtx and writes target genes (sparse TSV)

import os
import re
import time

t0 = time.time()
raw_dir = r"00_data\scRNA_cache\GSE176078_raw\Wu_etal_2021_BRCA_scRNASeq"
out_tsv = r"00_data\scRNA_cache\GSE176078_target16_counts.tsv"

target_genes = [
    "PPARG", "CXCL8", "IL1B", "FFAR2",
    "ADIPOQ", "PLIN1", "FABP4", "LEP",
    "EPCAM", "KRT19", "PTPRC", "CD68", "S100A8", "CSF3R", "VWF", "COL1A1",
]

with open(os.path.join(raw_dir, "count_matrix_genes.tsv"), "r", encoding="utf-8") as fg:
    genes = [line.strip() for line in fg]

with open(os.path.join(raw_dir, "count_matrix_barcodes.tsv"), "r", encoding="utf-8") as fb:
    barcodes = [line.strip() for line in fb]

gene_to_idx = {}
for idx1, g in enumerate(genes, start=1):
    if g in target_genes:
        gene_to_idx[str(idx1).encode("ascii")] = g

print(f"Loaded {len(genes)} genes and {len(barcodes)} barcodes; matched {len(gene_to_idx)}/{len(target_genes)} targets")
assert len(gene_to_idx) == len(target_genes), "Missing target genes!"

alt = b"|".join(sorted(gene_to_idx.keys(), key=lambda x: int(x)))
pattern = re.compile(rb"(?m)^(" + alt + rb")[ \t]+(\d+)[ \t]+(\d+)[ \t]*\r?$")

mtx_path = os.path.join(raw_dir, "count_matrix_sparse.mtx")
n_written = 0
chunk_size = 32 * 1024 * 1024

with open(mtx_path, "rb") as fm, open(out_tsv, "w", encoding="utf-8", buffering=1024 * 1024) as fo:
    fo.write("gene\tbarcode\tcount\n")
    hdr1 = fm.readline().decode("ascii", errors="replace").strip()
    hdr2 = fm.readline().decode("ascii", errors="replace").strip()
    print(f"MTX header: {hdr1} | dims: {hdr2}")
    parts = hdr2.split()
    assert int(parts[0]) == len(genes) and int(parts[1]) == len(barcodes), "Dimension mismatch!"

    buf = b""
    bytes_read = 0
    while True:
        chunk = fm.read(chunk_size)
        if not chunk:
            if buf:
                for m in pattern.finditer(buf):
                    g = gene_to_idx[m.group(1)]
                    bc = barcodes[int(m.group(2)) - 1]
                    val = m.group(3).decode("ascii")
                    fo.write(f"{g}\t{bc}\t{val}\n")
                    n_written += 1
            break
        bytes_read += len(chunk)
        data = buf + chunk
        nl = data.rfind(b"\n")
        if nl == -1:
            buf = data
            continue
        block = data[:nl + 1]
        buf = data[nl + 1:]
        out_lines = []
        for m in pattern.finditer(block):
            g = gene_to_idx[m.group(1)]
            bc = barcodes[int(m.group(2)) - 1]
            val = m.group(3).decode("ascii")
            out_lines.append(f"{g}\t{bc}\t{val}\n")
        if out_lines:
            fo.writelines(out_lines)
            n_written += len(out_lines)

dt = time.time() - t0
print(f"Done in {dt:.2f}s! Extracted {n_written} non-zero target entries -> {out_tsv}")
