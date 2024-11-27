#' @title delete bin on HiC matrix
#'
#' @description From DNA position of a deletion, the function remove the bins concern by the deletion
#'
#' @details
#' When a deletion start or end is positioned in a middle of a bin (i.e. not at the boundary), 3 solutions are possible:
#' -remove the bin from the matrix (limits = "max"),
#' -keep the bin in the matrix (limits = "min"),
#' -choose the nearest bin boundary (limits = "nearest").
#'
#'
#' @param matrix DNA input as DNAstring object (e.g chromosome sequence).
#' @param bin.width bin width of the matrix.
#' @param start one or more positions of deletion start in bp.
#' @param stop one or more positions of deletion stop/end in bp.
#' @param limits "nearest" to removed bins according to the nearest bin boundaries (default), max to removed bins involved in deletion, min to removed bins fully involved in deletion (see example and details).
#'
#' @return DNAstring
#'
#' @export
#'
#' @examples
#' matrix = sapply(1:5, function(x){rep(x, 5)}) ; matrix
#' bin.width = 10e3
#' start = 24e3
#' stop = 44e3
#'
#' #round deletion to range 20e3 and 40e3
#' delMAT(matrix, bin.width, start, stop, limits = "nearest")
#'
#' #round deletion to range 20e3 and 50e3
#' delMAT(matrix, bin.width, start, stop, limits = "max")
#'
#' #round deletion to range 30e3 and 40e3
#' delMAT(matrix, bin.width, start, stop, limits = "min")
#'
delMAT <- function(matrix, bin.width, start, stop, limits = "nearest") {

  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  bin.start.del = sapply(start, function(x) {
    ifelse(x %/% bin.width == x / bin.width, x / bin.width + 1,
           if (limits == "nearest") {round(x / bin.width) + 1} else if
           (limits == "max") {x %/% bin.width + 1} else if
           (limits == "min") {x %/% bin.width + 2})
  })

  bin.stop.del = sapply(stop, function(x) {
    ifelse(x %/% bin.width == x / bin.width, x / bin.width,
           if (limits == "nearest") {round(x / bin.width)} else if
           (limits == "max") {x %/% bin.width + 1} else if
           (limits == "min") {x %/% bin.width})
  })

  bin_to_removed = sapply(1:length(bin.start.del), function(i) {
    if (bin.start.del[i] <= bin.stop.del[i]) {bin.start.del[i]:bin.stop.del[i]}}
    ) %>% unlist

  bin_to_keep = (1:nrow(matrix))[!(1:nrow(matrix)) %in% bin_to_removed]

  return(matrix[bin_to_keep,bin_to_keep])

}

