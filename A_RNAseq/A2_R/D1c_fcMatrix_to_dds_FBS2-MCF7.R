#This runs DESEQ2 from featureCounts matrix ran on the command line version of subread
#It does the following (in sequential order):
#(1) first reads in the featureCounts txt file, filters for proteinCodingGenes only, saves as .rds
#(2) creates dds object from DESeqDataSetFromMatrix() and saves size-dispersion normalized mtrx as csv
#(3) runs DESeq() and saves paired DE results (l2fc, padj etc)
#Anney Ziqing Ye
#Last updated 7/30/2026
#Absolute parent directory paths have been anonymized (e.g. replaced with “/project_rootpath/”) for portability

library(DESeq2)

curr_proj_path <- "/home/username/curr_directory"
setwd(curr_proj_path) #set current working directory
source("helpers.R") #this path for source is relative to the working directory

#Make sure to add "/" at the end of outfilebasepath string
outfilebasepath <- "/project_rootpath/D1_DESeq2_R/"
codingonly_name <- "LY_FBS2-MCF7_fcCountMatrx_proteinCodingOnly.rds" #Enter desired output file name here. this is just a row-filtered unnormalized featureCounts matrix
fname_normcounts_human <- "LY_FBS2-MCF7_counts_DDSnormalized.csv" #Enter desired output file name here
ddsrds_name_human <- "LY_FBS2-MCF7_allConds_dds.rds" #This is mainly for use in other R scripts

# Read featureCounts table
fc_filepath <- "/project_rootpath/C1_featureCounts_cmd/LY_FBS2_MCF7_counts.txt"
fc <- read.delim(fc_filepath, comment.char = "#", check.names = FALSE)
gtf_file_path <- "/project_rootpath/Z_refGenome_files/gencode.v43.primary_assembly.annotation.gtf" #reference gtf file path 
#To ensure consistency, use the same exact gtf file as the one used in featureCounts

###################
# Part 1: read FeatureCounts table
countdata <- fc[, 7:ncol(fc)]  # Keep only count columns (first 6 columns are annotation)
rownames(countdata) <- fc$Geneid  # Use gene IDs as row names
# Optional: clean sample names if columns contain full BAM paths
colnames(countdata) <- basename(colnames(countdata)) #basename() automatically parses paths and gets filename
colnames(countdata) <- vapply(strsplit(colnames(countdata), "_", fixed = TRUE), \(p) p[2], character(1))
#Edit the line each time since the naming conventions may differ run-to-run
print(colnames(countdata))

# Make sure counts are integers #This is better done with assert statements
#print(is.integer(countdata)) #this'll return false cuz countdata is still dataframe
countdata_mtx <- as.matrix(countdata)
print(is.integer(countdata_mtx)) #should be true if the countdata slices only contained integers
stopifnot(typeof(countdata_mtx) == "integer",
          !anyNA(countdata_mtx),
          all(countdata_mtx >= 0)) #pre-requisite for DESeq2

host_counts  <- countdata_mtx

#discard rRNA and mtRNAs. Only include protein coding genes for DESeq2 analysis
#Make sure that the helpers.R file containing the following function is in the current directory
#and that it has been loaded via source("R/helpers.R")
fc_filtered <- filter_counts_mrna_no_rrna_no_mt(host_counts, gtf_file_path, return_extra = TRUE)

saveRDS(fc_filtered$counts, paste0(outfilebasepath, codingonly_name)) #Save a copy of the filtered FC table as RDS object


#########################
#Part 2. Run DESeq2 and generate dds object with normalized counts
countMatrix <- fc_filtered$counts #human counts

#Part 2.1. Enter sample metadata manually based on exact sample order in the fc header
treatment <- factor(
  c( "DMSO", "DMSO", "DMSO", "DMSO", "MA2", "MA2", "MA2", "MA2",
     "pIC", "pIC", "pIC", "pIC"),
  levels = c("DMSO", "MA2", "pIC") #DMSO equals to mock
)

colData <- data.frame(
  treatment = treatment,
  row.names = colnames(countMatrix)
)

#PART 2.2: Create DDS object, normalize based on size factors and save normalized counts
dds_human <- DESeqDataSetFromMatrix(
  countData = countMatrix, 
  colData = colData,
  design = ~ treatment
)
dds_human <- estimateSizeFactors(dds_human)
write.csv(counts(dds_human, normalized=TRUE), file = paste0(outfilebasepath, fname_normcounts_human))
host_sf <- sizeFactors(dds_human)
saveRDS(dds_human, paste0(outfilebasepath, ddsrds_name_human))


#PART 2.3: Differential gene expression analysis in DESeq 2 and save l2fc+padj results

write_deseq_results <- function(res, filename) { #helper function to write into csv file, res=results object
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  # Move gene_id to first column
  df <- df[, c("gene_id", setdiff(names(df), "gene_id"))]
  write.csv(df, file = filename, row.names = FALSE)
}

#Create files for each pair of treatment-to-mock comparison
dds <- DESeq(dds_human)

#2.3.1. DE MA2 vs mock (DMSO):
#results(dds, contrast = c("condition", "A", "M")) #general syntax comparing cond A against cond M
res1 <- results(dds, contrast = c("treatment", "MA2", "DMSO"))
write_deseq_results(res1, paste0(outfilebasepath, "LY_FBS2-MCF7_MA2_vs_mock.csv"))

#2.3.2.  DE diABZI vs mock (DMSO):
res2 <- results(dds, contrast = c("treatment", "pIC", "DMSO"))
write_deseq_results(res2, paste0(outfilebasepath, "LY_FBS2-MCF7_pIC_vs_mock.csv"))





