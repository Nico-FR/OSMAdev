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
#' @param DI logical. If TRUE, the function calculates the DI score. Default is FALSE.
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
#'   \item \strong{DI (Directionality Index / LFC Bias):} Computes the \eqn{\log_2} Fold Change of the directional
#'   bias between the MT and WT conditions for the specific mutated bin. It evaluates the shift in interactions
#'   towards the upstream (\eqn{U}) or downstream (\eqn{D}) regions.
#'   \deqn{\text{LFC}_{\text{bias}} = \frac{(D^{MT} - U^{MT}) - (D^{WT} - U^{WT})}{\log(2)}}
#'   Where \eqn{D} and \eqn{U} represent the mean interactions within a defined distance window from the mutated bin.
#' }
#'
#' Additionally, when multiple mutations occur within the exact same genomic coordinates (e.g., replicates or overlapping segments),
#' the function automatically appends logical columns (e.g., \code{max_corr_HFF}, \code{min_SIC_HFF}) to flag the mutation
#' exhibiting the strongest structural impact (lowest correlation or highest SIC) for downstream filtering.
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

analyseOrcaPredictions = function(predictions.dir, metadataWT, metadataMT, matrix.gz = FALSE, SIC = TRUE, corr = FALSE, DI = FALSE){

  ##############################
  # testing parameters
  #metadataWT = read.table("~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/metadataWT.tsv", header = TRUE, sep = "\t"); metadataMT = read.table("~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/metadataMT.tsv", header = TRUE, sep = "\t")
  #predictions.dir = "~/mnt/genome3D/Nicolas/inSilMut/chr1/step1/Predictions/" ; matrix.gz = TRUE ; SIC = TRUE ; corr = FALSE ; DI = FALSE
  ##############################

  if (unique(metadataMT$chr) %>% length != 1){
    stop("metadataMT must contain only one chromosome")
  }

  if (!SIC && !corr && !DI){
    stop("At least one of the parameters SIC, corr or DI must be TRUE")
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

  # Function 3: DI (LFC DI)
  fonction_DI <- function(WT.mat, MT.mat, start_mut, stop_mut, start_mat, bin.width, distanceBin) {
    if (DI) {

      # compute position of the mutated bin
      bin_mutated <- floor(
        ((start_mut + stop_mut) / 2 - start_mat) / bin.width
      ) + 1

      # Define Upstream and Downstream limits
      start_U <- bin_mutated - distanceBin
      end_U   <- bin_mutated - 1

      start_D <- bin_mutated + 1
      end_D   <- bin_mutated + distanceBin

      # Calculate mean U (Upstream / Left interactions)
      U_WT <- mean(WT.mat[start_U:end_U, bin_mutated], na.rm = TRUE)
      U_MT <- mean(MT.mat[start_U:end_U, bin_mutated], na.rm = TRUE)

      # Calculate mean D (Downstream / Right interactions)
      D_WT <- mean(WT.mat[bin_mutated, start_D:end_D], na.rm = TRUE)
      D_MT <- mean(MT.mat[bin_mutated, start_D:end_D], na.rm = TRUE)

      # Calculate Bias (Downstream - Upstream) = ln(Downstream / Upstream)
      Bias_WT <- D_WT - U_WT
      Bias_MT <- D_MT - U_MT

      # Calculate LFC Bias
      LFC_bias <- (Bias_MT - Bias_WT) / log(2)

      return(LFC_bias)
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
      WT.mat <- data.table::fread(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_hff", ext), header = FALSE, sep = "\t") %>% as.matrix()

      # ====================================================
      # corr OPTIMIZATION: Pre-calculate upper.tri once per WT matrix
      if (corr) {
        mask <- upper.tri(WT.mat)
        WT_up_tri <- WT.mat[mask]
      }

      # DI OPTIMIZATION: Pre-calculation of the DI parameters for this WT matrix
      if (DI) {
        bin.width <- metadataWT$scale[WT_idx] / 250
        start_mat <- metadataWT$start.mat[WT_idx]
        DI_distance <- WT_window - start_mat
        distanceBin <- floor(DI_distance / bin.width)
      }
      # ====================================================

      results_hff <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

        #get the ID of the MT matrix
        MT_row <- metadataMT.lst[[i]][i2, ]
        MT.ID <- MT_row$ID

        # Read the MT matrix using data.table::fread
        MT.mat <- data.table::fread(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_hff", ext), header = FALSE, sep = "\t") %>% as.matrix()

        list(
          corr_HFF = fonction_corr(WT_up_tri, MT.mat, mask),
          SIC_HFF  = fonction_SIC(WT.mat, MT.mat),
          DI_HFF   = fonction_DI(WT.mat, MT.mat, MT_row$start.mut, MT_row$stop.mut, start_mat, bin.width, distanceBin)
        )
      }) %>% dplyr::bind_rows()

      #select the results based on the parameters
      results_hff <- results_hff %>%
        dplyr::select(c(if(corr){1}, if(SIC){2}, if(DI){3}))

      # Combine the results into the metadataMT list
      metadataMT.lst[[i]] <- cbind(metadataMT.lst[[i]], results_hff)
    }



  ##############################
  # if ESC prediction
  ##############################
  if (metadataWT$model_ESC[WT_idx] == 1){

    # Read the WT matrix using data.table::fread (ULTRA FAST)
    WT.mat <- data.table::fread(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_esc", ext), header = FALSE, sep = "\t") %>% as.matrix()

    # ====================================================
    # corr OPTIMIZATION: Pre-calculate upper.tri once per WT matrix
    if (corr) {
      mask <- upper.tri(WT.mat)
      WT_up_tri <- WT.mat[mask]
    }

    # DI OPTIMIZATION: Pre-calculation of the DI parameters for this WT matrix
    if (DI) {
      bin.width <- metadataWT$scale[WT_idx] / 250
      start_mat <- metadataWT$start.mat[WT_idx]
      DI_distance <- WT_window - start_mat
      distanceBin <- floor(DI_distance / bin.width)
    }
    # ====================================================

    results_esc <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

      #get the ID of the MT matrix
      MT_row <- metadataMT.lst[[i]][i2, ]
      MT.ID <- MT_row$ID

      # Read the MT matrix using data.table::fread
      MT.mat <- data.table::fread(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_esc", ext), header = FALSE, sep = "\t") %>% as.matrix()

      list(
        corr_ESC = fonction_corr(WT_up_tri, MT.mat, mask),
        SIC_ESC  = fonction_SIC(WT.mat, MT.mat),
        DI_ESC   = fonction_DI(WT.mat, MT.mat, MT_row$start.mut, MT_row$stop.mut, start_mat, bin.width, distanceBin)
      )
    }) %>% dplyr::bind_rows()

    #select the results based on the parameters
    results_esc <- results_esc %>%
      dplyr::select(c(if(corr){1}, if(SIC){2}, if(DI){3}))

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
  }

  return(mutation_scores)
}



