#' @title Mutated metadata dataframe generation
#'
#' @description
#' The main goal of this function is to create MT metatdata to run Orca 1Mb predictions for all the 1Mb WT sequences.
#' For each 1Mb WT sequences, this function return the adjacent mutated windows that will be subject to mutations.
#'
#' @param metadataWT <dataFrame> : a dataframe generate by createMetadaWT
#' @param mutated.width <numeric> : size of the adjacent mutated window
#' @param rep <numeric> : the number of repetition for each mutations.
#'
#' @return dataframe
#'
#' @details
#' This function return a dataframe with metadatas relate to Mutated sequences used to run Orca predictions :
#' Column description:
#' * chr <numeric> : the number of the chromosome
#' * start.window <numeric> : start of the window that will contain the mutations
#' * stop.window <numeric> : stop of the window that will contain the mutations
#' * start.mat <numeric> : the first position of the 1Mb sequence
#' * stop.mat <numeric> : the last position of the 1Mb sequence
#' * start.mut <numeric> : start of mutated window
#' * stop.mut <numeric> : stop of the mutated window
#' * ID <character> : the name of the mutated type sequence "MT_x_y_z.fa" with: x the number of the related WildType sequence, y an incremental number for each different mutated window, z an incremental number for each repetition
#' * model = 1000000 : ORCA parameter cf ORCA documentation
#' * scale = 1000000 : ORCA parameter cf ORCA documentation
#' * wpos  500000 : ORCA parameter cf ORCA documentation
#' * mpos  500000 : ORCA parameter cf ORCA documentation
#' * model_HFF <0 or 1> : ORCA parameter. predict HFF matrices if 1
#' * model_ESC <0 or 1> : ORCA parameter. predict ESC matrices if 1
#'
#' @importFrom dplyr bind_rows arrange
#' @importFrom tibble tibble
#' @export


metadataMT = function(metadataWT, mutated.width, rep = 1){

  ####################################
  #testing parameters
  #metadataWT = metadat_WT ; mutated.width=4000 ; rep =1
  ####################################

  metadataMT = tibble::tibble()

  for(r in 1:rep){ # make the following processes as many time as rep value
    for (i in 1:nrow(metadataWT) ){ #for each line of the WildType metadata we make all the Mutated metadata line related.
      tmpTibble = tibble::tibble(
        chr = metadataWT$chr[i],
        start.window = metadataWT$start.window[i],
        stop.window = metadataWT$stop.window[i],
        start.mat = metadataWT$start.mat[i],
        stop.mat = metadataWT$stop.mat[i],
        start.mut = seq(metadataWT$start.window[i], metadataWT$stop.window[i] - mutated.width + 1, mutated.width),
        stop.mut = start.mut + mutated.width - 1,
        ID = paste0("MT_", i, "_", 1:length(start.mut), "_", r),
        model = 1e6,
        scale = 1e6,
        wpos = 0.5e6,
        mpos = 0.5e6,
        model_HFF = metadataWT$model_HFF[i],
        model_ESC = metadataWT$model_ESC[i]
      )
      # we merge all temporary tibble
      metadataMT = dplyr::bind_rows(metadataMT, after = tmpTibble)
    }

  }
  return(metadataMT %>% dplyr::arrange(start.mut))
}
