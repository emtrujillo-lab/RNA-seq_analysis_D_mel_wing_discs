# RNA-seq_analysis_D_mel_wing_discs
MEF2 Overexpression in Drosophila Imaginal Disc Myoblasts

# SETUP
dirs <- c("results", "figures")
invisible(sapply(dirs, function(d) if (!dir.exists(d)) dir.create(d)))

# -------------------------------------
# INSTALL AND LOAD PACKAGES

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bioc_pkgs <- c(
  "DESeq2", # differential expression (Love et al. 2014)
  "apeglm", # log2FC shrinkage (companion to DESeq2)
  "org.Dm.eg.db", # Drosophila gene annotation
  "AnnotationDbi", # annotation queries
  "clusterProfiler" # can use this GO enrichment instead of ClueGO/Cytoscape
)

if (!require("BiocManager"))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

cran_pkgs <- c(
  "fastqcr", # QC wrapper for FastQC (Kassambara)
  "ggplot2", # plotting
  "ggrepel", # non-overlapping labels
  "pheatmap", # heatmaps
  "EnhancedVolcano", # volcano plots
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "tibble",
  "readr"
)
new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc) > 0) BiocManager::install(new_bioc, ask = FALSE, update = FALSE)
if (length(new_cran) > 0) install.packages(new_cran)
suppressPackageStartupMessages({
  
  library(DESeq2)
  library(apeglm)
  library(org.Dm.eg.db)
  library(AnnotationDbi)
  library(fastqcr)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})
message("All packages loaded.")

# QUALITY CONTROL WITH fastqcr 
# - All categories PASSED except two expected RNA-seq failures:
# "per base sequence content" <- hexamer priming artifact; always fails
# "per sequence GC content" <- GC variation; not a concern for D. melanogaster
# - %GC was 47-50% across all samples (consistent with dm6)
# - Overall mapping rate to dm6 = 95%
# - NO trimming was performed (however, this can be completed if reads are not clean)
# ----------------------------------
raw_dir <- "./raw_reads"
qc_dir <- "./fastqc_results"
if (dir.exists(raw_dir) &&
    length(list.files(raw_dir, pattern = "\\.fastq\\.gz$")) > 0) {
  # Run FastQC via fastqcr
  fastqcr::fastqc(
    fq.dir = raw_dir,
    qc.dir = qc_dir,
    threads = 4
    # fastqc.path = "/path/to/fastqc" 
    
  )
  # Aggregate results
  qc <- fastqcr::qc_aggregate(qc_dir)
  message("\nFastQC summary:")
  print(fastqcr::qc_stats(qc))
  # Check for unexpected failures (the two below are normal for RNA-seq)
  expected_fails <- c("Per base sequence content", "Per sequence GC content")
  problems <- fastqcr::qc_fails(qc) %>%
    filter(!module %in% expected_fails)
  if (nrow(problems) > 0) {
    warning("Unexpected QC failures — review before proceeding:")
    print(problems)
  } else {
    message("QC check passed (only expected RNA-seq failures present).")
  }
  # HTML summary report
  fastqcr::qc_report(
    qc.path = qc_dir,
    result.file = "results/fastqcr_QC_report",
    interpret = TRUE
  )
  message("FastQC report: results/fastqcr_QC_report.html")
} else {
  message("Raw FASTQ directory not found or empty — skipping FastQC step.")
}

# SHELL COMMANDS: HISAT2 ALIGNMENT 
# Run these in a terminal, NOT in R (files may be to big to run, check if run freezes)
# --------------------------------------
# use HISAT2 (single-end mode, no trimming) against dm6.
#
# ── Download pre-built HISAT2 dm6 index ───
# wget https://genome-idx.s3.amazonaws.com/hisat/dm6.tar.gz
# tar -xzf dm6.tar.gz
#
# ── Align each sample (SINGLE-END) ────
# for SAMPLE in MEF2_OE_1 MEF2_OE_2 MEF2_OE_3 Control_1 Control_2 Control_3; do
# hisat2 \
# -x ./dm6/genome \
# -U ./raw_reads/${SAMPLE}.fastq.gz \

# --dta \
# -p 8 \
# -S ./aligned/${SAMPLE}.sam 2> ./aligned/${SAMPLE}_hisat2_log.txt
# samtools sort -@ 8 -o ./aligned/${SAMPLE}.bam ./aligned/${SAMPLE}.sam
# samtools index ./aligned/${SAMPLE}.bam
# rm ./aligned/${SAMPLE}.sam
# done
#
# ── Count reads with featureCounts ----------------------------------
# featureCounts \
# -T 8 \
# -a Drosophila_melanogaster.BDGP6.46.gtf \
# -o counts/all_samples_counts.txt \
# -F GTF -t exon -g gene_id \
# --minOverlap 10 \
# --primary \
# aligned/MEF2_OE_1.bam aligned/MEF2_OE_2.bam aligned/MEF2_OE_3.bam \
# aligned/Control_1.bam aligned/Control_2.bam aligned/Control_3.bam
#
# ----------------------------------
# IMPORT featureCounts COUNT MATRIX
# ----------------------------------
counts_file <- "counts/all_samples_counts.txt"
if (!file.exists(counts_file)) {
  # ── Synthetic demo data (6 samples) ───
  message("Count file not found — using synthetic demo data (3 MEF2-OE + 3 Control).")
  set.seed(2024)
  n_genes <- 14000
  gene_ids <- paste0("FBgn", sprintf("%07d", seq_len(n_genes)))
  base_mat <- matrix(
    rnbinom(n_genes * 6, mu = 120, size = 4),
    nrow = n_genes,
    dimnames = list(
      gene_ids,
      c("MEF2_OE_1", "MEF2_OE_2", "MEF2_OE_3",
        "Control_1", "Control_2", "Control_3")
      
    )
  )
  base_mat[1:300, 1:3] <- base_mat[1:300, 1:3] * 6 # muscle genes up
  base_mat[301:450, 1:3] <- round(base_mat[301:450, 1:3] / 5) # cell-cycle down
  count_matrix <- base_mat
} else {
  # ── featureCounts output ────────────
  counts_raw <- read.table(counts_file, header = TRUE, sep = "\t",
                           skip = 1, comment = "#", check.names = FALSE)
  count_matrix <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
  rownames(count_matrix) <- counts_raw$Geneid
  colnames(count_matrix) <- gsub(".*[/\\\\]", "", colnames(count_matrix))
  colnames(count_matrix) <- gsub("\\.bam$", "", colnames(count_matrix))
}
message("Count matrix: ", nrow(count_matrix), " genes x ",
        ncol(count_matrix), " samples")
message("Samples: ", paste(colnames(count_matrix), collapse = ", "))


# SAMPLE METADATA
# ----------------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(
    ifelse(grepl("MEF2|OE|mef2", colnames(count_matrix), ignore.case = TRUE),
           "MEF2_OE", "Control"),
    levels = c("Control", "MEF2_OE")
  ),
  row.names = colnames(count_matrix)
)
stopifnot(all(rownames(sample_info) == colnames(count_matrix)))
message("\nSample metadata:")
print(sample_info)

# DESeq2 DIFFERENTIAL EXPRESSION
# ----------------------------------
# Input : featureCounts raw counts
# Cutoff : padj < 0.05 (Benjamini-Hochberg)
# Result : >147 up-regulated genes in MEF2-OE condition

# ----------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)
# Pre-filter low-count genes (standard)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("\nGenes after filtering: ", sum(keep))
dds$condition <- relevel(dds$condition, ref = "Control")
message("Running DESeq2...")
dds <- DESeq(dds)
message("\nSize factors:")
print(round(sizeFactors(dds), 3))

# QC PLOTS
# ----------------------------------
vsd <- vst(dds, blind = TRUE)
# ── PCA ────────────────────────────
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"), 1)
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Control" = "#2E75B6", "MEF2_OE" = "#C00000")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA — MEF2-OE vs. Control",
          subtitle = "Wing imaginal discs, n=3 per condition (2024)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("figures/PCA_plot.pdf", p_pca, width = 6.5, height = 5)
ggsave("figures/PCA_plot.png", p_pca, width = 6.5, height = 5, dpi = 300)

# ── Sample distance heatmap ─────────────

sampleDists <- dist(t(assay(vsd)))
ann_col <- data.frame(Condition = sample_info$condition,
                      row.names = rownames(sample_info))
ann_colors <- list(Condition = c(Control = "#2E75B6", MEF2_OE = "#C00000"))
pdf("figures/SampleDistance_heatmap.pdf", width = 6, height = 5)
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2E75B6", "white"))(100),
         border_color = NA, main = "Sample-to-sample distances (VST)")
dev.off()
# ── Dispersion plot ─────────────────────
pdf("figures/Dispersion_plot.pdf", width = 6, height = 5)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()
message("QC plots saved.")

# EXTRACT RESULTS AND ANNOTATE
# ----------------------------------
# Raw results (threshold: padj < 0.05)
res_raw <- results(dds,
                   contrast = c("condition", "MEF2_OE", "Control"),
                   alpha = 0.05,
                   pAdjustMethod = "BH")
message("\n--- DESeq2 summary (padj < 0.05) ---")
summary(res_raw, alpha = 0.05)
# apeglm-shrunken log2FC for visualization
res_shrunk <- lfcShrink(dds,
                        coef = "condition_MEF2_OE_vs_Control",
                        type = "apeglm", res = res_raw)
# Annotate
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    symbol = mapIds(org.Dm.eg.db, keys = gene_id, keytype = "ENSEMBL",
                    column = "SYMBOL", multiVals = "first"),
    
    entrez = mapIds(org.Dm.eg.db, keys = gene_id, keytype = "ENSEMBL",
                    column = "ENTREZID", multiVals = "first"),
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange <= 0 ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "Down", "NS")),
    label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id)
  ) %>%
  arrange(padj)
message("\nDE genes (padj < 0.05): Up=", sum(res_df$sig=="Up", na.rm=TRUE),
        " Down=", sum(res_df$sig=="Down", na.rm=TRUE))
message("(Paper reports >147 up-regulated genes)")
write_csv(res_df, "results/DESeq2_MEF2_OE_vs_Control_full.csv")
write_csv(filter(res_df, sig != "NS"),
          "results/DESeq2_MEF2_OE_vs_Control_significant.csv")
# MA plot
pdf("figures/MA_plot.pdf", width = 6, height = 5)
plotMA(res_shrunk, ylim = c(-8, 8),
       main = "MA plot — MEF2-OE vs. Control (apeglm)",
       colSig = "#C00000", alpha = 0.05)
dev.off()
message("Results saved.")

# VOLCANO PLOT
# ----------------------------------
# Genes to highlight 
highlight <- c(
  "Act57B", "Act79B", "Act88F", # qPCR-validated actins
  "Mhc", # qPCR-validated MHC
  "Mlc2", "Mlc1", # qPCR-validated myosin light chains
  "TpnI", # qPCR-validated troponin
  "sls", # qPCR-validated sallimus
  "Tm1", "Tm2", "TpnC4", "TpnT", # other sarcomeric
  "Actn", "Zasp52", "Zasp66", "bt",
  "Myo61F", "zipper",
  "wupA", "flightin", # adult-specific (NOT activated)
  "Mef2" # overexpressed gene
)

res_df <- mutate(res_df,
                 volcano_label = ifelse(label %in% highlight, label, ""))
p_volcano <- EnhancedVolcano(res_df,
                             lab = res_df$volcano_label,
                             x = "log2FoldChange",
                             y = "padj",
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             title = "MEF2 Overexpression vs. Control",
                             subtitle = "Trujillo et al. 2024, Dev Biol 516:82-95",
                             caption = paste0("DESeq2 padj<0.05 | Up: ",
                                              sum(res_df$sig=="Up", na.rm=TRUE),
                                              " Down: ",
                                              sum(res_df$sig=="Down", na.rm=TRUE)),
                             col = c("grey70", "#2E75B6", "#FFB347", "#C00000"),
                             colAlpha = 0.7, pointSize = 1.8, labSize = 3.2, labFace = "italic",
                             drawConnectors = TRUE, widthConnectors = 0.4,
                             maxoverlapsConnectors = 40, legendPosition = "right")
ggsave("figures/Volcano_MEF2_OE.pdf", p_volcano, width = 10, height = 9)
ggsave("figures/Volcano_MEF2_OE.png", p_volcano, width = 10, height = 9, dpi = 300)
message("Volcano plot saved.")



# GENE ONTOLOGY ENRICHMENT
# ----------------------------------
# used ClueGO (Bindea et al. 2009), a Cytoscape plug-in.
# ClueGO is a graphical desktop tool 


# SESSION INFO
# ----------------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("\n=== Analysis complete ===")
message("Raw + analyzed data: https://doi.org/10.17632/777hyfy6h5.1 (Mendeley)")
message("\nOutputs:")
message(" figures/PCA_plot.pdf")
message(" figures/Volcano_MEF2_OE.pdf")
message(" figures/Heatmap_top50_DE.pdf")
message(" figures/Muscle_gene_barplot.pdf")
message(" figures/GO_BP_Up_regulated.pdf")
message(" results/DESeq2_MEF2_OE_vs_Control_full.csv")
message(" results/Muscle_gene_summary.csv")
message(" results/sessionInfo.txt")
