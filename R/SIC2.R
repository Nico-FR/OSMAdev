#' @title Structural impact score (bin wise)  between 2 matrices
#'
#' @description Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between each bin of 2 matrices.
#' Fold change is computed for each bin of the matrices provided
#'
#' @details The function takes as input 2 matrices (mutated and wildtype) and computes the average absolute log fold change between each bins (return one value for the entire matrices).
#'
#'
#' @param mutated.mat,wildtype.mat mutated (query) or wildtype (control, subject) matrices (observed or observed/expected matrices) as `dgCMatrix` or `matrix` object for only one chromosome.
#'
#' @importFrom magrittr %>%
#'
#' @return integer
#'
#' @export
#'
#' @examples
#' # to do
#'
SIC <- function(mutated.mat, wildtype.mat) {

  . <- NULL

  #sanity check
  if(!inherits(mutated.mat, c("Matrix", "matrix"))) {
    stop("mutated matrix is not a matrix or Matrix object")}
  if(!inherits(wildtype.mat, c("Matrix", "matrix"))) {
    stop("wildtype matrix is not a matrix or Matrix object")}
  if (nrow(wildtype.mat) != nrow(mutated.mat)){
    warning("matrices do not have the same size!")
  }
  if (ncol(wildtype.mat) != ncol(mutated.mat)){
    warning("matrices do not have the same size!")
  }

  # make sure matrices are symmetric
  mutated.mat <- as.matrix(mutated.mat)
  wildtype.mat <- as.matrix(wildtype.mat)

  if (!isSymmetric(mutated.mat)) {
    mutated.mat[lower.tri(mutated.mat)] <- t(mutated.mat)[lower.tri(mutated.mat)]
  }
  if (!isSymmetric(wildtype.mat)) {
    wildtype.mat[lower.tri(wildtype.mat)] <- t(wildtype.mat)[lower.tri(wildtype.mat)]
  }

  #remove 0
  mutated.mat[mutated.mat == 0] <- NA
  wildtype.mat[wildtype.mat == 0] <- NA

  # Calculate SIC
  fold_change <- mutated.mat / wildtype.mat
  SIC <- fold_change |> log() |> abs() |> Matrix::rowMeans(na.rm = TRUE)

  return(SIC)

}






















