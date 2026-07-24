#' @title Decrease the resolution of an interaction matrix
#'
#' @description `matZoomify` decreases the resolution of a Hi-C interaction matrix by binning.
#' For example, converting a matrix from a 50kb resolution to a 100kb resolution.
#'
#' @details
#' The input matrix is first converted to a sparse `CsparseMatrix` and restricted to its upper triangular part (including the diagonal).
#' Then, bins are aggregated by summing the interaction counts according to the ratio of `new_bin_width` to `old_bin_width`.
#' If the original number of bins is not a multiple of the reduction factor, a final smaller bin is created at the end to contain the remaining bins.
#' The resulting matrix is returned as an upper triangular `dgCMatrix` (sparse Matrix).
#'
#' @param matrix `dgCMatrix` or `matrix` object representing intra-chromosomal interactions.
#' @param old_bin_width Bin width of the input matrix in base pairs (i.e. original resolution).
#' @param new_bin_width Target bin width in base pairs (i.e. new resolution). Must be a multiple of `old_bin_width`.
#' @param verbose Logical. Whether or not to print messages. Default = `TRUE`.
#'
#' @return An upper triangular `dgCMatrix` object with the decreased resolution.
#'
#' @importFrom Matrix triu t sparseMatrix
#' @importFrom methods as
#' @export
#'
#' @examples
#' # Decrease resolution of HCT116 chromosome 19 matrix from 50kb to 100kb
#' mat_100kb <- matZoomify(mat_HCT116_chr19_50kb,
#'                         old_bin_width = 50e3,
#'                         new_bin_width = 100e3)
#'
matZoomify <- function(matrix, old_bin_width, new_bin_width, verbose = TRUE) {

  # 1. Sanity checks
  if (!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")
  }

  if (old_bin_width <= 0) {
    stop("old_bin_width must be positive")
  }

  if (new_bin_width <= 0) {
    stop("new_bin_width must be positive")
  }

  if (new_bin_width %% old_bin_width != 0) {
    stop("new_bin_width must be a multiple of old_bin_width")
  }

  factor <- as.integer(new_bin_width / old_bin_width)

  if (factor == 1) {
    if (verbose) message("new_bin_width is equal to old_bin_width, returning original matrix as dgCMatrix.")
    return(methods::as(methods::as(Matrix::triu(matrix), "CsparseMatrix"), "generalMatrix"))
  }

  n_bins <- nrow(matrix)

  if (verbose) {
    message(paste0("Reducing resolution from ", old_bin_width, " bp to ", new_bin_width,
                   " bp (reduction factor = ", factor, ")."))
    message(paste0("Original matrix dimensions: ", n_bins, "x", n_bins))
  }

  # 2. Convert to sparse and restrict to upper triangular part
  m_sparse <- methods::as(matrix, "CsparseMatrix")
  m_upper <- Matrix::triu(m_sparse)

  # 3. Bin aggregation using S matrix
  n_new <- ceiling(n_bins / factor)
  group_idx <- rep(1:n_new, each = factor, length.out = n_bins)

  S <- Matrix::sparseMatrix(i = group_idx, j = 1:n_bins, x = 1)
  res <- S %*% m_upper %*% Matrix::t(S)
  res_upper <- Matrix::triu(res)
  res_out <- methods::as(res_upper, "generalMatrix")

  if (verbose) {
    message(paste0("New matrix dimensions: ", nrow(res_out), "x", ncol(res_out)))
  }

  return(res_out)
}
