#' @title Log2 mean fold change between 2 matrices
#'
#' @description Fold change is computed for each bin of the matrices provided and the absolute log2 mean fold change is returned.
#'
#' @param mat1,mat2 input matrices as `dgCMatrix` or `matrix` object for only one chromosome.
#'
#' @return integer
#'
#' @export
#'
#' @examples
#' # to do
#'
SIC <- function(mat1, mat2) {

  #sanity check
  if(!inherits(mat1, c("Matrix", "matrix"))) {
    stop("mat1 is not a matrix or Matrix object")}

  if(!inherits(mat2, c("Matrix", "matrix"))) {
    stop("mat2 is not a matrix or Matrix object")}

  if (nrow(mat2) != nrow(mat1) | ncol(mat2) != ncol(mat1)) {
    warning("matrices do not have the same size!")
  }

  # get upper triangle + diag values as vectors
  vect1 <- mat1[upper.tri(mat1, diag = TRUE)]
  vect2 <- mat2[upper.tri(mat2, diag = TRUE)]

  #remove 0
  vect1[vect1 == 0] <- NA
  vect2[vect2 == 0] <- NA

  # Calculate SIC
  fold_change <- vect1 / vect2
  SIC <- fold_change |> log2 () |> abs() |> mean(na.rm = TRUE)

  return(SIC)

}
