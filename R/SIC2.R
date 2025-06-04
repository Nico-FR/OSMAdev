#' @title Structural impact score between 2 matrices
#'
#' @description Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between 2 matrices.
#' Fold change is computed for each bin of the matrices provided
#'
#' @param mutated.mat,wildtype.mat mutated (query) or wildtype (control, subject) matrices as `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Bin width of the matrix in base pair.
#' @param vp.start Start of the view point in base pair.
#' @param vp.stop Stop/end of the view point in base pair. Default is NULL to only use the bin where vp.start is located (ie: vp.stop = vp.start + bin.width).
#' @param start,stop Area in bp of the chromosome to compute impact score. Default is NULL to use the entire chromosome (i.e. entire matrix).
#' @param verbose if TRUE, show information when view point does not match a bin annotations.
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
SIC <- function(mutated.mat, wildtype.mat, bin.width, vp.start, vp.stop = NULL, start = NULL, stop = NULL, verbose = TRUE) {

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



  if (!isSymmetric(mutated.mat)) {
    mutated.mat[lower.tri(mutated.mat)] <- t(mutated.mat)[lower.tri(mutated.mat)]
  }
  if (!isSymmetric(wildtype.mat)) {
    wildtype.mat[lower.tri(wildtype.mat)] <- t(wildtype.mat)[lower.tri(wildtype.mat)]
  }

  # Calculate SIC
  diff = mutated.mat - wildtype.mat
  SIC <- diff %>% abs %>% Matrix::rowMeans(., na.rm = TRUE)

  return(SIC)

}






















