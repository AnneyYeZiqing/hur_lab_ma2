

## 1. Reusable function for importing GENCODE v43 GTF, keep protein_coding mRNAs,
##    explicitly drop rRNA and MT-RNA


filter_counts_mrna_no_rrna_no_mt <- function(count_mat,
                                             gtf_file,
                                             keep_gene_types = "protein_coding",
                                             exclude_gene_types = c("rRNA"),
                                             exclude_mt_by_name = TRUE,
                                             mt_prefix = "^MT-",
                                             return_extra = FALSE,
                                             verbose = TRUE) {
  stopifnot(is.matrix(count_mat) || is.data.frame(count_mat))
  count_mat <- as.matrix(count_mat)
  stopifnot(!is.null(rownames(count_mat)))
  stopifnot(file.exists(gtf_file))
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it with install.packages('data.table').")
  }
  
  # Read GTF (skip comment lines)
  gtf_df <- data.table::fread(
    cmd = paste("grep -v '^#' ", shQuote(gtf_file)),
    sep = "\t",
    header = FALSE
  )
  colnames(gtf_df) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )
  
  # Keep only gene rows 
  genes_gtf <- gtf_df[gtf_df$feature == "gene", , drop = FALSE]
  
  # Safe attribute extractor for GENCODE-style "key \"value\";"
  extract_attr <- function(attr_vec, key) {
    sapply(strsplit(attr_vec, ";"), function(fields) {
      fields <- trimws(fields)
      hit <- fields[startsWith(fields, paste0(key, " "))]
      if (length(hit) == 0) return(NA_character_)
      sub(paste0('^', key, ' "?([^"]*)"?'), '\\1', hit[1])
    }, USE.NAMES = FALSE)
  }
  
  gene_annot <- data.frame(
    gene_id   = extract_attr(genes_gtf$attribute, "gene_id"),
    gene_type = extract_attr(genes_gtf$attribute, "gene_type"),
    gene_name = extract_attr(genes_gtf$attribute, "gene_name"),
    stringsAsFactors = FALSE
  )
  
  # Map annotation to count matrix rows
  idx <- match(rownames(count_mat), gene_annot$gene_id)
  gene_type <- gene_annot$gene_type[idx]
  gene_name <- gene_annot$gene_name[idx]
  
  # Build keep mask
  keep <- !is.na(gene_type) & gene_type %in% keep_gene_types
  if (!is.null(exclude_gene_types) && length(exclude_gene_types) > 0) {
    keep <- keep & !(gene_type %in% exclude_gene_types)
  }
  if (exclude_mt_by_name) {
    keep <- keep & !is.na(gene_name) & !grepl(mt_prefix, gene_name)
  }
  
  ## ---- Added: compute and print per-sample mtRNA / rRNA / unknown percentages BEFORE dropping ----
  if (verbose) {
    n_samps <- ncol(count_mat)
    sample_names <- colnames(count_mat)
    if (is.null(sample_names)) sample_names <- paste0("Sample", seq_len(n_samps))
    
    # per-sample totals (columns)
    total_per_sample <- colSums(count_mat, na.rm = TRUE)
    
    # masks for rows
    mt_mask <- !is.na(gene_name) & grepl(mt_prefix, gene_name)
    rrna_mask <- !is.na(gene_type) & (gene_type %in% exclude_gene_types)
    unknown_mask <- is.na(idx)  # rows in count_mat that didn't match a gene_id in GTF
    
    # per-sample read sums for each mask (handle no TRUE rows safely)
    if (any(mt_mask)) {
      mt_per_sample <- colSums(count_mat[mt_mask, , drop = FALSE], na.rm = TRUE)
    } else {
      mt_per_sample <- rep(0, n_samps)
    }
    if (any(rrna_mask)) {
      rrna_per_sample <- colSums(count_mat[rrna_mask, , drop = FALSE], na.rm = TRUE)
    } else {
      rrna_per_sample <- rep(0, n_samps)
    }
    if (any(unknown_mask)) {
      unknown_per_sample <- colSums(count_mat[unknown_mask, , drop = FALSE], na.rm = TRUE)
    } else {
      unknown_per_sample <- rep(0, n_samps)
    }
    
    # compute percentages (safe for zero totals)
    pct_mt <- ifelse(total_per_sample > 0, 100 * mt_per_sample / total_per_sample, NA_real_)
    pct_rrna <- ifelse(total_per_sample > 0, 100 * rrna_per_sample / total_per_sample, NA_real_)
    pct_unknown <- ifelse(total_per_sample > 0, 100 * unknown_per_sample / total_per_sample, NA_real_)
    
    # Print a concise per-sample summary line
    cat("Pre-filter per-sample read summary (total, mtRNA, rRNA, unknown):\n")
    for (i in seq_len(n_samps)) {
      cat(sprintf(
        "%s: total=%d; mtRNA=%d (%.2f%%); rRNA=%d (%.2f%%); unknown=%d (%.2f%%)\n",
        sample_names[i],
        total_per_sample[i],
        mt_per_sample[i], ifelse(is.na(pct_mt[i]), 0, pct_mt[i]),
        rrna_per_sample[i], ifelse(is.na(pct_rrna[i]), 0, pct_rrna[i]),
        unknown_per_sample[i], ifelse(is.na(pct_unknown[i]), 0, pct_unknown[i])
      ))
    }
    cat("\n")
  }
  ## ---- end added block ----
  
  if (verbose) {
    if (any(is.na(idx))) {
      message("Unmatched gene IDs (not found in GTF gene rows): ", sum(is.na(idx)))
    }
    message("Keeping genes: ", sum(keep), " / ", length(keep))
    message("Kept gene_type breakdown:")
    print(table(gene_type[keep], useNA = "ifany"))
    if (exclude_mt_by_name) {
      message("Mito genes remaining (should be FALSE only):")
      print(table(grepl(mt_prefix, gene_name[keep]), useNA = "ifany"))
    }
  }
  
  filtered <- count_mat[keep, , drop = FALSE]
  
  if (return_extra) {
    return(list(
      counts = filtered,
      keep = keep,
      gene_type = gene_type,
      gene_name = gene_name,
      gene_annot = gene_annot
    ))
  } else {
    return(filtered)
  }
}


# Usage examples: 
# If you want extra info (keep mask + annotations)
#out <- filter_counts_mrna_no_rrna_no_mt(count_mat, gtf_file, return_extra = TRUE)
#count_mat_filtered <- out$counts

