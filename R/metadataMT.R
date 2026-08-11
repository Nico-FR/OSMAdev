#' @title Mutated metadata dataframe generation
#'
#' @description
#' The main goal of this function is to create MT metadata to create fasta files for each mutation.
#' For each 1Mb WT sequences, this function return all the adjacent mutation positions on the mutated window.
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
#' * start.pred <numeric> : the first position of the 1Mb sequence
#' * stop.pred <numeric> : the last position of the 1Mb sequence
#' * start.mut <numeric> : start of mutated window
#' * stop.mut <numeric> : stop of the mutated window
#' * ID <character> : the name of the mutated type sequence "MT_x_y_z.fa" with: x the number of the related WildType sequence, y an incremental number for each different mutated window, z an incremental number for each repetition
#'
#' @importFrom dplyr bind_rows arrange
#' @importFrom tibble tibble
#' @export


metadataMT = function(metadataWT, mutated.width, rep = 1){

  ####################################
  #testing parameters
  #metadataWT = metadat_WT ; mutated.width=4000 ; rep =1
  ####################################

  metadataMT_list = list()

  for(r in 1:rep){ # make the following processes as many time as rep value
    metadataMT_list[[r]] = lapply(1:nrow(metadataWT), function(i) {
      base_df <- tibble::tibble(
        chr = metadataWT$chr[i],
        start.window = metadataWT$start.window[i],
        stop.window = metadataWT$stop.window[i],
        start.pred = metadataWT$start.pred[i],
        stop.pred = metadataWT$stop.pred[i],
        start.mut = seq(metadataWT$start.window[i], metadataWT$stop.window[i] - mutated.width + 1, mutated.width),
        stop.mut = start.mut + mutated.width - 1,
        ID = paste0("MT_", i, "_", 1:length(start.mut), "_", r)
      )

      other_cols <- setdiff(names(metadataWT), c("chr", "start.window", "stop.window", "start.pred", "stop.pred", "ID", "A", "T", "G", "C"))
      other_df <- metadataWT[i, other_cols, drop = FALSE]
      other_df_recycled <- other_df[rep(1, nrow(base_df)), , drop = FALSE]
      row.names(other_df_recycled) <- NULL

      tibble::as_tibble(cbind(base_df, other_df_recycled))
    }) %>% dplyr::bind_rows()
  }
  metadataMT = dplyr::bind_rows(metadataMT_list)
  return(metadataMT %>% dplyr::arrange(start.mut))
}
