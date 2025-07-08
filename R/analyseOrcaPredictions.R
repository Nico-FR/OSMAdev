#' @title Calculate scores to measure mutation effects
#'
#' @description
#' This function calculates scores (see details) between WT matrices and their corresponding MT matrices.
#' It reads the prediction files from a specified directory, processes them based on the provided metadatas, and computes the scores for each mutant matrix against its corresponding wild-type matrix.
#' The function supports both HFF and ESC models, and can handle gzipped matrix files.
#' The function return the metadataMT with the scores added as new columns.
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
#' The function calculates the following scores:
#' - **SIC**: The mean absolute log fold change between the MT and WT matrices.
#' - **corr**: The Pearson correlation coefficient between the MT and WT matrices.
#' - **DI**: Currently not implemented, returns NA.
#'
#' @return data.frame
#'
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom progress progress_bar
#' @importFrom dplyr bind_rows select group_by mutate ungroup
#' @importFrom utils read.table
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

  # Fonction 1: corr
  fonction_corr <- function(WT, MT) {
    if (corr) {
      return(stats::cor(MT[upper.tri(MT)], WT[upper.tri(WT)], use = "pairwise.complete.obs", method = "pearson"))
    } else {return(NA)}
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


  # Fonction 2: SIC
  fonction_SIC <- function(WT, MT) {
    if (SIC) {
      diff = MT - WT
      return(diff %>% abs %>% mean(na.rm = TRUE))
    } else {return(NA)}
  }

  # Fonction 3: DI
  fonction_DI <- function(WT, MT) {
    if (DI) {
      return(NA)
    } else {return(NA)}
  }

  #split metadataMT by start.window
  metadataMT.lst <- split(metadataMT, metadataMT$start.window)

  #initiate progress bar
  pb <- progress::progress_bar$new(total = length(metadataMT.lst),
                                   format = "  calculating [:bar] :percent in :elapsed",
                                   clear = FALSE, width= 60)
  pb$tick(0)

  #loop for each start.mat in metadataMT.lst (i.e each MT matrix)
  for (i in 1:length(metadataMT.lst)){

    #get the corresponding WT matrix
    WT.ID = metadataWT$ID[metadataWT$start.window == as.numeric(names(metadataMT.lst[i]))]

    ###############################
    # if HFF prediction
    ###############################
    if (metadataWT$model_HFF[i] == 1){

      #read the WT matrix
      if (matrix.gz){
        WT.mat = read.table(gzfile(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_hff.tsv.gz")), header = FALSE, sep = "\t") %>% as.matrix()
      } else {
        WT.mat = read.table(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_hff.tsv"), header = FALSE, sep = "\t") %>% as.matrix()
      }

      results_hff <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

        #get the ID of the MT matrix
        MT.ID = metadataMT.lst[[i]]$ID[i2]

        #read the MT matrix
        if (matrix.gz){
          MT.mat = read.table(gzfile(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_hff.tsv.gz")), header = FALSE, sep = "\t") %>% as.matrix()
        } else {
          MT.mat = read.table(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_hff.tsv"), header = FALSE, sep = "\t") %>% as.matrix()
        }

        data.frame(
          corr_HFF = fonction_corr(WT.mat, MT.mat),
          SIC_HFF = fonction_SIC(WT.mat, MT.mat),
          DI_HFF = fonction_DI(WT.mat, MT.mat)
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
  if (metadataWT$model_ESC[i] == 1){

    #read the WT matrix
    if (matrix.gz){
      WT.mat = read.table(gzfile(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_esc.tsv.gz")), header = FALSE, sep = "\t") %>% as.matrix()
    } else {
      WT.mat = read.table(paste0(predictions.dir, WT.ID, "_predictions_1000000_1M_esc.tsv"), header = FALSE, sep = "\t") %>% as.matrix()
    }

    results_esc <- lapply(1:nrow(metadataMT.lst[[i]]), function(i2) {

      #get the ID of the MT matrix
      MT.ID = metadataMT.lst[[i]]$ID[i2]

      #read the MT matrix
      if (matrix.gz){
        MT.mat = read.table(gzfile(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_esc.tsv.gz")), header = FALSE, sep = "\t") %>% as.matrix()
      } else {
        MT.mat = read.table(paste0(predictions.dir, MT.ID, "_predictions_1000000_1M_esc.tsv"), header = FALSE, sep = "\t") %>% as.matrix()
      }

      data.frame(
        corr_ESC = fonction_corr(WT.mat, MT.mat),
        SIC_ESC = fonction_SIC(WT.mat, MT.mat),
        DI_ESC = fonction_DI(WT.mat, MT.mat)
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
  mutation_scores = do.call(rbind, metadataMT.lst)
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



