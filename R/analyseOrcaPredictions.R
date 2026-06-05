#' @title Calculate scores to measure mutation effects
#'
#' @description
#' This function calculates scores to measure the effect of mutations on chromatin structure based on WT matrices (\eqn{mat_{WT} = log{\frac{Observed}{Expected}}}) and their corresponding MT matrices (\eqn{mat_{MT} = log{\frac{Observed}{Expected}}}).
#' It reads the prediction files from a specified directory, processes them based on the provided metadatas, and computes the scores (see details).
#'
#' @param predictions.dir character. The directory containing the prediction files.
#' @param metadataWT data.frame. The metadata for the wild-type matrices.
#' @param metadataMT data.frame. The metadata for the mutant matrices.
#' @param matrix.gz logical. If TRUE, the function expects gzipped matrix files. Default is FALSE.
#' @param SIC logical. If TRUE, the function calculates the SIC score. Default is TRUE.
#' @param corr logical. If TRUE, the function calculates the correlation score. Default is FALSE.
#' @param localSIC logical. If TRUE, the function calculates the local SIC score (lSIC). Default is FALSE.
#'
#' @details
#' The scores are calculated for each mutant (MT) matrix against its corresponding wild-type (WT) matrix,
#' based on the provided metadata. The function handles both HFF and ESC models seamlessly.
#'
#' \itemize{
#'   \item \strong{SIC (Structural Impact Score):} Measures the mean absolute logarithmic fold change
#'   between the MT and WT matrices. Since the input matrices are assumed to be already \eqn{\log(\text{Observed} / \text{Expected})} transformed,
#'   the fold change is computed as the simple difference between the two matrices.
#'   \deqn{\text{SIC} = \frac{1}{N} \sum_{(i,j) \in \Omega} \left| M^{MT}_{i,j} - M^{WT}_{i,j} \right|}
#'   Where \eqn{N} is the total number of valid bins in the upper triangle (\eqn{\Omega}).
#'
#'   \item \strong{corr (Pearson Correlation):} Calculates the Pearson correlation coefficient between the upper
#'   triangular values of the MT and WT matrices, excluding missing values (pairwise complete observations).
#'   \deqn{r = \frac{\sum (M^{MT}_{i,j} - \bar{M}^{MT})(M^{WT}_{i,j} - \bar{M}^{WT})}{\sqrt{\sum (M^{MT}_{i,j} - \bar{M}^{MT})^2 \sum (M^{WT}_{i,j} - \bar{M}^{WT})^2}}}
#'
#'   \item \strong{lSIC (local Structural Impact Score):} Measures the mean absolute logarithmic fold change
#'   between the MT and WT matrices for the mutated bin and its adjacent neighbors (+/- 1 bin).
#'   It evaluates the impact within a defined distance window, excluding the diagonal.
#'   \deqn{\text{lSIC} = \frac{1}{N} \sum_{(i,j) \in \Omega_{local}} \left| M^{MT}_{i,j} - M^{WT}_{i,j} \right|}
#' }
#'
#' Additionally, when multiple mutations occur within the exact same genomic coordinates (i.e., replicates),
#' the function automatically appends logical columns (e.g., \code{max_corr_HFF}, \code{min_SIC_HFF}, \code{min_lSIC_HFF}) to flag the mutation
#' exhibiting the strongest (or weakest) structural impact for downstream filtering.
#'
#' @return data.frame
#'
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom progress progress_bar
#' @importFrom dplyr bind_rows select group_by mutate ungroup
#' @importFrom data.table fread
#'
#' @export

analyseOrcaPredictions = function(predictions.dir, metadataWT, metadataMT, matrix.gz = FALSE, SIC = TRUE, corr = FALSE, localSIC = FALSE){

  ##############################
  # testing parameters
  #metadataWT = read.table("~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/metadataWT.tsv", header = TRUE, sep = "\t"); metadataMT = read.table("~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/metadataMT.tsv", header = TRUE, sep = "\t")
  #predictions.dir = "~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/Predictions/" ; matrix.gz = TRUE ; SIC = TRUE ; corr = FALSE ; localSIC = FALSE
  ##############################

  if (unique(metadataMT$chr) %>% length != 1){
    stop("metadataMT must contain only one chromosome")
  }

  if (!SIC && !corr && !localSIC){
    stop("At least one of the parameters SIC, corr or localSIC must be TRUE")
  }

  #check prediction files
  ID.path = list.files(predictions.dir, pattern = "predictions_1000000") %>%
    gsub("_predictions_1000000_1M_hff.tsv.gz", "", .) %>%
    gsub("_predictions_1000000_1M_esc.tsv.gz", "", .)

  ID.metadatas = c(metadataWT$ID, metadataMT$ID)

  if (!all(ID.path %in% ID.metadatas)) {
    warning("Some prediction files have no corresponding IDs in metadataWT and metadataMT.")
  }

  if (!all(ID.metadatas %in% ID.path)) {
    warning("Some IDs in metadataWT or metadataMT do not have corresponding prediction files in the predictions directory.")
  }

  # Fonction 1: corr
  fonction_corr <- function(WT_up_tri, MT, mask) {
    if (corr) {
      return(stats::cor(MT[mask], WT_up_tri, use = "pairwise.complete.obs", method = "pearson"))
    } else {return(NA)}
  }

  # Function 2: SIC
  fonction_SIC <- function(WT, MT) {
    if (SIC) {
      return(mean(abs(MT - WT), na.rm = TRUE))
    } else {return(NA)}
  }

  # Function 3: lSIC (local SIC)
  fonction_lSIC <- function(WT.mat, MT.mat, start_mut, stop_mut, start_mat, bin.width, distanceBin) {
    if (localSIC) {

      # compute position of the mutated bin
      bin_mutated <- floor(
        ((start_mut + stop_mut) / 2 - start_mat) / bin.width
      ) + 1

      # Define the range of bins (+/- 1)
      bins_target <- (bin_mutated - 1):(bin_mutated + 1)
      # Ensure bins are within matrix bounds
      bins_target <- bins_target[bins_target >= 1 & bins_target <= ncol(WT.mat)]

      # Calculate absolute differences
      diff_mat <- abs(MT.mat - WT.mat)

      # For each target bin, collect its interactions within distanceBin
      values <- c()
      for (b in bins_target) {
        start_idx <- max(1, b - distanceBin)
        end_idx <- min(ncol(WT.mat), b + distanceBin)

        # Get indices in this row, excluding the diagonal (row == col)
        indices <- start_idx:end_idx
        indices <- indices[indices != b]

        values <- c(values, diff_mat[b, indices])
      }

      return(mean(values, na.rm = TRUE))
    } else {
      return(NA)
    }
  }

  #split metadataMT by start.window
  metadataMT.lst <- split(metadataMT, metadataMT$start.window)

  #initiate progress bar
  pb <- progress::progress_bar$new(total = length(metadataMT.lst),
                                   format = "  calculating [:bar] :percent in :elapsed",
                                   clear = FALSE, width= 60)
  pb$tick(0)

  #Dynamically managed extension for fread
  ext <- if(matrix.gz) ".tsv.gz" else ".tsv"

  #loop for each start.mat in metadataMT.lst (i.e each MT matrix)
  for (i in 1:length(metadataMT.lst)){

    #get the corresponding WT matrix
    WT_window <- as.numeric(names(metadataMT.lst[i]))
    WT_idx <- which(metadataWT$start.window == WT_window)[1]
    WT.ID <- metadataWT$ID[WT_idx]

    ###############################
    # if HFF prediction
    ###############################
    if (metadataWT$model_HFF[WT_idx] == 1){

      # Read the WT matrix using data.table::fread (ULTRA FAST)
      WT.mat <- data.table::fread(file.path(predictions.dir, paste0(WT.ID, "_predictions_1000000_1M_hff", ext)), header = FALSE, sep = "\t") %>% as.matrix()

      # ====================================================
      # corr OPTIMIZATION: Pre-calculate upper.tri once per WT matrix
      if (corr) {
        mask <- upper.tri(WT.mat)
        WT_up_tri <- WT.mat[mask]
      }

      # lSIC OPTIMIZATION: Pre-calculation of the lSIC parameters for this WT matrix
      if (localSIC) {
        bin.width <- metadataWT$scale[WT_idx] / 250
        start_mat <- metadataWT$start.mat[WT_idx]
        lSIC_distance <- WT_window - start_mat
        distanceBin <- floor(lSIC_distance / bin.width)
      }
      # ====================================================

      results_hff <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

        #get the ID of the MT matrix
        MT_row <- metadataMT.lst[[i]][i2, ]
        MT.ID <- MT_row$ID

        # Read the MT matrix using data.table::fread
        MT.mat <- data.table::fread(file.path(predictions.dir, paste0(MT.ID, "_predictions_1000000_1M_hff", ext)), header = FALSE, sep = "\t") %>% as.matrix()

        list(
          corr_HFF = fonction_corr(WT_up_tri, MT.mat, mask),
          SIC_HFF  = fonction_SIC(WT.mat, MT.mat),
          lSIC_HFF = fonction_lSIC(WT.mat, MT.mat, MT_row$start.mut, MT_row$stop.mut, start_mat, bin.width, distanceBin)
        )
      }) %>% dplyr::bind_rows()

      #select the results based on the parameters
      results_hff <- results_hff %>%
        dplyr::select(c(if(corr){1}, if(SIC){2}, if(localSIC){3}))

      # Combine the results into the metadataMT list
      metadataMT.lst[[i]] <- cbind(metadataMT.lst[[i]], results_hff)
    }



  ##############################
  # if ESC prediction
  ##############################
  if (metadataWT$model_ESC[WT_idx] == 1){

    # Read the WT matrix using data.table::fread (ULTRA FAST)
    WT.mat <- data.table::fread(file.path(predictions.dir, paste0(WT.ID, "_predictions_1000000_1M_esc", ext)), header = FALSE, sep = "\t") %>% as.matrix()

    # ====================================================
    # corr OPTIMIZATION: Pre-calculate upper.tri once per WT matrix
    if (corr) {
      mask <- upper.tri(WT.mat)
      WT_up_tri <- WT.mat[mask]
    }

    # lSIC OPTIMIZATION: Pre-calculation of the lSIC parameters for this WT matrix
    if (localSIC) {
      bin.width <- metadataWT$scale[WT_idx] / 250
      start_mat <- metadataWT$start.mat[WT_idx]
      lSIC_distance <- WT_window - start_mat
      distanceBin <- floor(lSIC_distance / bin.width)
    }
    # ====================================================

    results_esc <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

      #get the ID of the MT matrix
      MT_row <- metadataMT.lst[[i]][i2, ]
      MT.ID <- MT_row$ID

      # Read the MT matrix using data.table::fread
      MT.mat <- data.table::fread(file.path(predictions.dir, paste0(MT.ID, "_predictions_1000000_1M_esc", ext)), header = FALSE, sep = "\t") %>% as.matrix()

      list(
        corr_ESC = fonction_corr(WT_up_tri, MT.mat, mask),
        SIC_ESC  = fonction_SIC(WT.mat, MT.mat),
        lSIC_ESC = fonction_lSIC(WT.mat, MT.mat, MT_row$start.mut, MT_row$stop.mut, start_mat, bin.width, distanceBin)
      )
    }) %>% dplyr::bind_rows()

    #select the results based on the parameters
    results_esc <- results_esc %>%
      dplyr::select(c(if(corr){1}, if(SIC){2}, if(localSIC){3}))

    # Combine the results into the metadataMT list
    metadataMT.lst[[i]] <- cbind(metadataMT.lst[[i]], results_esc)

  }
    pb$tick()

  }

  # Combine the results into a single data frame
  mutation_scores = dplyr::bind_rows(metadataMT.lst)
  row.names(mutation_scores) <- NULL

  # best scores analysis in the case of multiple mutations in the same region (i.e. rep >= 2 : multiple mutations with the same start.mut and stop.mut)
  if(any(duplicated(mutation_scores$start.mut))) {#if there is replicates (same region mutated more than ones)
    # add TRUE or FALSE to indicate the highest (SIC) or lowest (corr) score
    if ("corr_HFF" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(max_corr_HFF = corr_HFF == max(corr_HFF, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    if ("SIC_HFF" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(min_SIC_HFF = SIC_HFF == min(SIC_HFF, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    if ("corr_ESC" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(max_corr_ESC = corr_ESC == max(corr_ESC, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    if ("SIC_ESC" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(min_SIC_ESC = SIC_ESC == min(SIC_ESC, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    if ("lSIC_HFF" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(min_lSIC_HFF = lSIC_HFF == min(lSIC_HFF, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    if ("lSIC_ESC" %in% colnames(mutation_scores)) {
      mutation_scores <- mutation_scores %>%
        dplyr::group_by(start.mut, stop.mut) %>%
        dplyr::mutate(min_lSIC_ESC = lSIC_ESC == min(lSIC_ESC, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
  }

  return(mutation_scores)
}



