#' @title Structural impact score between 2 matrices
#'
#' @description Measurement described in Zhou 2022 use to compute the average absolute log fold change of interactions between 2 matrices.
#' Fold change can be compute:
#'
#' * from one bin to all bins in the matrix (using view point parameters),
#'
#' * to one bin to all bins in the matrix (using start/stop parameters).
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


  #matrix coordinate of limits
  if (is.null(start)) {
    l.start <- 1}
  if (is.null(stop)) {
    l.stop <- nrow(mutated.mat)}

  if (!is.null(start)) {
    l.start <- start %/% bin.width + 1}
  if (!is.null(stop)) {
    l.stop <- stop %/% bin.width}

  #matrix coordinate of view point
  if (is.null(vp.stop)) {vp.stop <- (vp.start %/% bin.width + 1) * bin.width} #vp.stop = (bin.start + 1 bin) * bin.width

  if ((vp.start %/% bin.width == vp.start / bin.width) & (vp.stop %/% bin.width == vp.stop / bin.width)) {
    vp.start <- vp.start / bin.width + 1 #+1 because bin 1 start to 0 and stop to bin.width
    vp.stop <- vp.stop / bin.width
  } else {
    vp.start <- vp.start %/% bin.width + 1 # +1 because matrix coordinates start to 1. So first nucleotide (i.e. vp.start = 1bp) is bin nb 1.
    vp.stop <- ifelse(vp.stop %/% bin.width <= vp.start, vp.start, vp.stop %/% bin.width)
    if (verbose == TRUE) {
      message(paste0("vp.start/vp.stop are not multiples of bin.width, round to ", (vp.start - 1) * bin.width,
                   " and ", vp.stop * bin.width, " (i.e. ", vp.stop + 1 - vp.start, " bins)."))
    }
  }

  #crop matrices to area
  MT <- as.matrix(mutated.mat[l.start:l.stop, l.start:l.stop])
  if (!isSymmetric(MT)) {
    MT[lower.tri(MT)] <- t(MT)[lower.tri(MT)]
  }
  WT <- as.matrix(wildtype.mat[l.start:l.stop, l.start:l.stop])
  if (!isSymmetric(WT)) {
    WT[lower.tri(WT)] <- t(WT)[lower.tri(WT)]
  }

  #crop matrices to vp
  MT <- MT[(vp.start - l.start + 1):(vp.stop - l.start + 1),]
  WT <- WT[(vp.start - l.start + 1):(vp.stop - l.start + 1),]

  # Calculate SIC
  diff = MT - WT
  SIC <- diff %>% abs %>% mean(., na.rm = TRUE)

  return(SIC)

}






















