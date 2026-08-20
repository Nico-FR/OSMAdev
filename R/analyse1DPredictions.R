#' @title Calculate fold change to measure mutation effects on 1D track predictions
#'
#' @description
#' This function calculates fold change scores to measure the effect of mutations on 1D genomic signal predictions
#' (such as RNA-seq expression levels, ChIP-seq histone/TF signal, CAGE, etc.).
#' It reads predicted 1D track values (BedGraph files) from a specified directory, processes them based on WT and MT metadatas,
#' and computes global and optional local fold changes across all predicted tracks.
#'
#' @param predictions.dir character. The directory containing the prediction files (BedGraph format) and the two metadata files.
#' @param metadataWT data.frame. The metadata for the wild-type sequences.
#' @param metadataMT data.frame. The metadata for the mutant sequences.
#' @param regions.gr GRanges. Optional GRanges object representing specific genomic regions of interest to calculate local fold change. Default is NULL.
#' @param bedgraph.gz logical. Default is FALSE. If TRUE, the function searches for gzipped files (.bedgraph.gz).
#' @param track_nums numeric. Optional vector of track numbers to analyze (see AlphagenomePredictionsMetadata.tsv). If specified, only tracks with track_num in track_nums are processed. Default is NULL.
#'
#' @return A data.table containing the mutation metadata with appended fold change scores (gFC_x, min_gFC_x) for all tracks.
#'
#' @importFrom progress progress_bar
#' @importFrom data.table fread fwrite setDT set as.data.table
#' @importFrom GenomicRanges start end seqnames
#'
#' @export
analyse1DPredictions = function(predictions.dir, metadataWT, metadataMT, regions.gr = NULL, bedgraph.gz = FALSE, track_nums = NULL) {

  # -------------------------------------------------------------------------
  # 1. INPUT VALIDATION & METADATA PARSING
  # -------------------------------------------------------------------------

  # Ensure all mutations target the same chromosome to avoid coordinate mismatches
  if (length(unique(metadataMT$chr)) != 1){
    stop("metadataMT must contain only one chromosome")
  }

  if (!dir.exists(predictions.dir)) {
    stop("predictions.dir does not exist: ", predictions.dir)
  }

  # Identify AlphaGenome metadata files generated during the prediction step
  meta_files <- list.files(predictions.dir, pattern = "AlphagenomePredictionsMetadata\\.tsv$", full.names = TRUE)
  if (length(meta_files) != 2) {
    stop("Expected exactly 2 files matching 'AlphagenomePredictionsMetadata.tsv' in ", predictions.dir, " but found ", length(meta_files))
  }

  # Load metadata files using data.table::fread for fast I/O
  meta1 <- data.table::fread(meta_files[1], header = TRUE, sep = "\t")
  meta2 <- data.table::fread(meta_files[2], header = TRUE, sep = "\t")

  # Verify structural integrity: both metadata files must share key tracking columns
  cols_to_check <- intersect(c("name", "track_num", "ontology_curie"), colnames(meta1))
  if (length(cols_to_check) == 0) {
    stop("Metadata files lack required track columns (name, track_num, or ontology_curie).")
  }

  # Note: The '..' prefix in data.table tells R to look for 'cols_to_check' in the
  # calling environment rather than treating it as a column name inside the table itself.
  m1_sub <- meta1[, ..cols_to_check]
  m2_sub <- meta2[, ..cols_to_check]

  if (!identical(m1_sub, m2_sub)) {
    stop("The two metadata files matching 'AlphagenomePredictionsMetadata.tsv' are not identical.")
  }

  suffix <- ifelse(isTRUE(bedgraph.gz), ".bedgraph.gz", ".bedgraph")

  # -------------------------------------------------------------------------
  # 2. TRACK FILTERING
  # -------------------------------------------------------------------------

  # If the user requested specific tracks, filter the metadata accordingly
  if (!is.null(track_nums)) {
    meta1 <- meta1[track_num %in% track_nums]
    meta2 <- meta2[track_num %in% track_nums]

    if (nrow(meta1) == 0) {
      stop("No tracks matched the specified track_nums: ", paste(track_nums, collapse = ", "))
    }
  }

  # Build unique identifiers for each track to map them to physical files
  track_ids <- paste0(meta1$ontology_curie, "_", meta1$track_num)
  message("Analysing ", length(track_ids), " predicted tracks: ", paste(track_ids, collapse = ", "))

  # -------------------------------------------------------------------------
  # 3. PRE-ALLOCATION & COORDINATE PRE-COMPUTATION
  # -------------------------------------------------------------------------

  # Convert the output table to data.table format
  results_df <- data.table::as.data.table(metadataMT)

  # Pre-allocate gFC and lFC columns for all tracks
  gfc_cols <- c()
  lfc_cols_by_track <- list()

  # for each tracks IDs
  for (j in seq_along(track_ids)) {
    track_num <- meta1$track_num[j]
    gfc_col <- paste0("gFC_", track_num)

    # Create gFC column (We use NA_real_ instead of NA to strictly enforce the 'numeric' type natively).
    data.table::set(results_df, j = gfc_col, value = NA_real_)
    gfc_cols <- c(gfc_cols, gfc_col)

    lfc_cols_by_track[[j]] <- character()
    if (!is.null(regions.gr) && length(regions.gr) > 0) {
      if (!inherits(regions.gr, "GRanges")) {
        stop("regions.gr must be a GRanges object")
      }

      # Generate column names for local fold change based on region names (default) or boundaries
      region_names <- names(regions.gr)
      lfc_cols_by_track[[j]] <- sapply(seq_along(regions.gr), function(k) {
        r_start <- GenomicRanges::start(regions.gr)[k]
        r_end <- GenomicRanges::end(regions.gr)[k]

        label <- if (!is.null(region_names) && !is.na(region_names[k]) && region_names[k] != "") {
          region_names[k]
        } else {
          paste0(r_start, "_", r_end)
        }
        paste0("lFC_", label, "_", track_num)
      })

      # Pre-allocate local Fold Change columns
      for (lfc_col in lfc_cols_by_track[[j]]) {
        data.table::set(results_df, j = lfc_col, value = NA_real_)
      }
    }
  }

  # Precompute relative coordinates for each region and each MT row to optimize speed
  region_rel_coords <- list()
  if (!is.null(regions.gr) && length(regions.gr) > 0) {
    for (k in seq_along(regions.gr)) {
      r_chr <- as.character(GenomicRanges::seqnames(regions.gr)[k])
      r_start <- GenomicRanges::start(regions.gr)[k]
      r_end <- GenomicRanges::end(regions.gr)[k]

      # Determine overlapping coordinates between the region and the predicted window
      overlap_start <- pmax(r_start, results_df$start.pred)
      overlap_end <- pmin(r_end, results_df$stop.pred)

      overlaps <- (results_df$chr == r_chr) & (overlap_start <= overlap_end)

      # Convert absolute genomic coordinates into relative indices for the 1D prediction array
      rel_starts <- ifelse(overlaps, overlap_start - results_df$start.pred + 1, NA_integer_)
      rel_ends <- ifelse(overlaps, overlap_end - results_df$start.pred + 1, NA_integer_)

      region_rel_coords[[k]] <- list(
        rel_start = rel_starts,
        rel_end = rel_ends,
        overlaps = overlaps
      )
    }
  }

  # -------------------------------------------------------------------------
  # 4. BATCH PROCESSING LOOP
  # -------------------------------------------------------------------------

  # Create a temporary column 'orig_idx' using '.I' (data.table's internal row number).
  # This guarantees we update the correct row in 'results_df', even after splitting the data.
  results_df[, orig_idx := .I]

  # Split the workload by start.window. This allows us to load the Wild Type (WT)
  # file only once per window, significantly reducing expensive disk I/O operations.
  metadataMT.lst <- split(results_df, results_df$start.window)

  pb <- progress::progress_bar$new(
    total = nrow(results_df),
    format = "  processing mutation :current/:total [:bar] :percent in :elapsed",
    clear = FALSE, width = 60
  )
  pb$tick(0)

  for (i in seq_along(metadataMT.lst)) {

    MT_sub <- metadataMT.lst[[i]]
    WT_window <- as.numeric(names(metadataMT.lst[i]))
    WT_idx <- which(metadataWT$start.window == WT_window)[1]

    if (is.na(WT_idx)) {
      warning("No matching WT window found for start.window: ", WT_window)
      pb$tick(nrow(MT_sub))
      next
    }

    WT.ID <- metadataWT$ID[WT_idx]

    # Pre-load all WT track arrays for the current window into memory.
    # Reading files is the slowest part of the loop, so we batch it here.
    WT_values_list <- vector("list", length(track_ids))
    for (j in seq_along(track_ids)) {
      track_id <- track_ids[j]
      WT_file <- file.path(predictions.dir, paste0(WT.ID, "_predictions_1048576_", track_id, suffix))
      if (file.exists(WT_file)) {
        # '[[1]]' instantly extracts the first column as a raw numeric vector
        WT_values_list[[j]] <- data.table::fread(WT_file, header = FALSE, sep = "\t")[[1]]
      }
    }

    # Iterate through every Mutant (MT) sequence tied to the current WT window
    for (r in seq_len(nrow(MT_sub))) {
      MT.ID <- MT_sub$ID[r]
      orig_r <- MT_sub$orig_idx[r] # The true row number in the master results_df table

      for (j in seq_along(track_ids)) {
        track_id <- track_ids[j]
        track_num <- meta1$track_num[j]
        gfc_col <- paste0("gFC_", track_num)

        WT_values <- WT_values_list[[j]]

        MT_file <- file.path(predictions.dir, paste0(MT.ID, "_predictions_1048576_", track_id, suffix))

        if (!file.exists(MT_file)) {
          warning("MT prediction file not found: ", MT_file)
          next
        }

        # Load the corresponding Mutant values array
        MT_values <- data.table::fread(MT_file, header = FALSE, sep = "\t")[[1]]

        # COMPUTE GLOBAL FC
        gFC_val <- mean(MT_values - WT_values, na.rm = TRUE)

        # IMPORTANT OPTIMIZATION: data.table::set does not create a copy of the dataframe
        # bypassing the usual memory bottleneck generated by results_df[[gfc_col]][orig_r] <- gFC_val.
        data.table::set(results_df, i = orig_r, j = gfc_col, value = gFC_val)

        # COMPUTE LOCAL FC (if specific GRanges were provided)
        if (!is.null(regions.gr) && length(regions.gr) > 0) {
          for (k in seq_along(regions.gr)) {
            lfc_col <- lfc_cols_by_track[[j]][k]
            coords <- region_rel_coords[[k]]

            # Check if the mutation physically falls within the region bounds
            if (coords$overlaps[orig_r]) {
              rel_start <- coords$rel_start[orig_r]
              rel_end <- coords$rel_end[orig_r]

              lFC_val <- mean(MT_values[rel_start:rel_end] - WT_values[rel_start:rel_end], na.rm = TRUE)
              data.table::set(results_df, i = orig_r, j = lfc_col, value = lFC_val)
            }
          }
        }
      }
      pb$tick() # Update progress bar after each MT sequence
    }
  }

  # Clean up: remove the temporary mapping index as we no longer need it
  results_df[, orig_idx := NULL]

  # -------------------------------------------------------------------------
  # 5. REPLICATE FLAGS (MINIMUM FC EVALUATION)
  # -------------------------------------------------------------------------

  if (any(duplicated(results_df$start.mut))) {

    # Gather all newly calculated Fold Change columns
    target_cols <- c(gfc_cols, unlist(lfc_cols_by_track))

    # Generate the min_ column names using same label logic
    min_cols <- paste0("min_", target_cols)

    # VECTORIZED GROUP-BY OPERATION:
    # Instead of looping over each track with dplyr, we evaluate all columns at once.
    # - '.SDcols = target_cols' tells data.table which columns to pass to the function.
    # - '.SD' (Subset of Data) represents these columns.
    # - 'by = .(start.mut, stop.mut)' groups the rows by mutation locus natively.
    # - 'lapply' runs the minimum check across all target columns simultaneously in C.
    results_df[, (min_cols) := lapply(.SD, function(x) {
      if (all(is.na(x))) return(NA)
      x == min(x, na.rm = TRUE)
    }), by = .(start.mut, stop.mut), .SDcols = target_cols]
  }

  return(results_df)
}
