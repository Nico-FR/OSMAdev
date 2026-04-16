#' @title Structural impact score (bin wise)  between 2 matrices
#'
#' @description Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between each bin of 2 matrices.
#' Fold change is computed for each bin of the matrices provided
#'
#' @details The function takes as input 2 matrices (mutated and wildtype) and computes the average absolute log fold change between each bins.
#'
#'
#' @param mat1,mat2 input matrices as `dgCMatrix` or `matrix` object for only one chromosome.
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
SIC <- function(mat1, mat2) {

  . <- NULL

  #sanity check
  if(!inherits(mat1, c("Matrix", "matrix"))) {
    stop("mutated matrix is not a matrix or Matrix object")}
  if(!inherits(mat2, c("Matrix", "matrix"))) {
    stop("wildtype matrix is not a matrix or Matrix object")}
  if (nrow(mat2) != nrow(mat1)){
    warning("matrices do not have the same size!")
  }
  if (ncol(mat2) != ncol(mat1)){
    warning("matrices do not have the same size!")
  }

  # make sure matrices are symmetric
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)

  if (!isSymmetric(mat1)) {
    mat1[lower.tri(mat1)] <- t(mat1)[lower.tri(mat1)]
  }
  if (!isSymmetric(mat2)) {
    mat2[lower.tri(mat2)] <- t(mat2)[lower.tri(mat2)]
  }

  #remove 0
  mat1[mat1 == 0] <- NA
  mat2[mat2 == 0] <- NA

  # Calculate SIC
  fold_change <- mat1 / mat2
  SIC <- fold_change |> log2() |> abs() |> Matrix::rowMeans(na.rm = TRUE)

  return(SIC)

}






















