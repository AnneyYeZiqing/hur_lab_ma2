###################################
# PCA from DESeq2-normalized DDS
#Absolute parent directory paths have been anonymized (e.g. replaced with “/project_rootpath/”) for portability

library(DESeq2)
library(ggplot2)
library(ggrepel)


basepath <- "/project_rootpath/D1_DESeq2_R/"
dds_rds_path <- "/project_rootpath/D1_DESeq2_R/LY_FBS10-HCC2429_allConds_dds.rds"

# ---- Addl parameters (edit paths as needed) ----
out_pdf        <- paste0(basepath, "PCA_vst_PC1_vs_PC2_LY-FBS10-HCC2429.pdf") #syntax for concat string
out_png        <- paste0(basepath, "PCA_vst_PC1_vs_PC2_LY-FBS10-HCC2429.png")
pca_coords_csv <- paste0(basepath, "PCA_vst_PC1_vs_PC2_LY-FBS10-HCC2429.csv")


if (!file.exists(dds_rds_path)) stop("dds RDS not found: ", dds_rds_path)
dds <- readRDS(dds_rds_path)
if (!inherits(dds, "DESeqDataSet")) stop("Loaded object is not a DESeqDataSet.")


message("Loaded DESeqDataSet with ", nrow(dds), " genes and ", ncol(dds), " samples.")

# ---- Ensure size factors (needed for vst normalization) ----
if (is.null(sizeFactors(dds))) {
  message("Size factors not found. Estimating size factors...")
  dds <- estimateSizeFactors(dds)
} else {
  message("Found existing size factors.")
}

# dds = your DESeqDataSet
# Make sure you have a sample-name column in colData:
# e.g., colData(dds)$sample == "WT_mock1", "WT_mock2", ...
stopifnot("treatment" %in% names(colData(dds)))

vsd <- vst(dds, blind = TRUE)

# ---- Extract transformed matrix (genes x samples) ----
mat_vst <- assay(vsd)
# Confirm orientation
if (!all(colnames(mat_vst) == colnames(dds))) {
  warning("Column names of transformed matrix and dds differ; proceeding but check sample alignment.")
}

# ---- PCA using ALL genes (no gene selection) ----
# prcomp expects samples x features -> transpose genes x samples -> samples x genes
mat_for_pca <- t(mat_vst)
message("Performing PCA on all genes (", ncol(mat_for_pca), " features). This may be slower for very large datasets.")

pca_res <- prcomp(mat_for_pca, center = TRUE, scale. = FALSE) # do not scale after VST

###################################
# #Alternatively: PCA using filtered genes:
# # --- PCA on top variable genes ---
# library(matrixStats)
# rv <- matrixStats::rowVars(assay(vsd))
# n_top <- 1000
# top <- order(rv, decreasing = TRUE)[seq_len(min(n_top, length(rv)))]
# 
# mat_for_pca <- t(assay(vsd)[top, , drop = FALSE])
# pca_res <- prcomp(mat_for_pca, center = TRUE, scale. = FALSE)

##############end of alternative block #######################

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# ---- Prepare PCA coordinates and metadata ----
pca_df <- as.data.frame(pca_res$x)  # rows are samples (rownames = sample IDs used in PCA)

cd <- as.data.frame(colData(vsd))

# Align metadata to PCA row order (CRITICAL)
cd <- cd[rownames(pca_df), , drop = FALSE]

# use colnames as sample labels
sample_code <- colnames(dds)
print(sample_code)

pca_df$sample <- sample_code

# Combine (now guaranteed aligned)
pca_df <- cbind(pca_df, cd)

# Save coordinates
write.csv(pca_df, pca_coords_csv, row.names = FALSE)
message("Saved PCA coordinates to: ", pca_coords_csv)

#################################
# Plot
# ---- Plot PC1 vs PC2 ----
pcx <- 1 #if want to plot PC3, 4, etc, change the numbers accordingly
pcy <- 2
xlab <- sprintf("PC%d (%.1f%%)", pcx, var_explained[pcx])
ylab <- sprintf("PC%d (%.1f%%)", pcy, var_explained[pcy])

# Build plot in a simple, robust way (no conditional inside aes)
if ("treatment" %in% colnames(pca_df)) { #"treatment" name corresponds to those specified in dds object during dds construction
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment))
} else {
  p <- ggplot(pca_df, aes(x = PC1, y = PC2))
}

p <- p +
  geom_point(size = 3) +
  labs(x = xlab, y = ylab, title = "PCA") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Add labels only if sample column exists
if ("sample" %in% colnames(pca_df)) {
  p <- p + ggrepel::geom_text_repel(aes(label = sample),
                                    size = 3, show.legend = FALSE)
}
print(sample_code)
colnames(pca_df)

ggsave(out_pdf, plot = p, width = 7, height = 6)
ggsave(out_png, plot = p, width = 7, height = 6, dpi = 300)
message("Saved PCA plot to: ", out_pdf, " and ", out_png)