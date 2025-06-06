#' @title Import Orca matrix into full chromosome matrix
#'
#' @description This function create an empty matrix of the full chromosome which is then filled with Orca matrix (i.e. 250x250 bins).
#'
#' @details From the matrices produced by orca_predict.genomepredict function.
#'
#' @param df_prediction.path The full path of the Orca matrix. Data frame without header which contains orca predictions, i.e. log(Obseved/Expected).
#' @param scale Size of the scale (i.e. window) of the predictions in base pair (i.e. 32Mb, 16Mb, 8Mb, 4Mb, 2Mb and 1Mb).
#' @param mpos The coordinate to zoom into for multiscale prediction. By default mpos is the center of the 32Mb sequence provided (i.e. start + 16Mb).
#' @param chromsize Chromosome size in base pair of the related chromosome.
#' @param output Character: Default is "OE" to return observed / expected counts. "Obs" to return observed counts.
#' @param df_normmats.path Optional. Background distance-based expected balanced contact matrices. Only needed to return observed counts.
#' @param sep separator between matrix field (df_prediction.path or df_normmats.path), default is tabulation.
#' @param model orca prediction model (32e6 or 1e6). Default is 32e6.
#'
#' @return A `Matrix` class object: upper triangular and sparse Matrix
#'
#' @importFrom Matrix triu sparseMatrix
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#'
#' @export
#'
#'
orca2matrix <- function(df_prediction.path, sep = "\t", mpos, scale, chromsize, output = "OE", df_normmats.path = NULL, model = 32e6) {

  #matrix specifications
  bin.width = scale / 250
  nbins = ifelse(chromsize %/% bin.width == chromsize / bin.width, chromsize %/% bin.width, chromsize %/% bin.width + 1)  #nb bins of the final matrix

  start.tmp = ifelse(mpos - scale / 2 < 0, 0, mpos - scale / 2) # theoretical start in bp (>= 0)

  #shift (- 1 bin) start.tmp if model != scale
  start = ifelse(scale == model, start.tmp,
                 start.tmp - bin.width) #(start.tmp %/% (bin.width * 2) + 1) * (bin.width * 2) - (bin.width * 2))  #start.tmp %/% (bin.width * 2) * (bin.width * 2)
  if(start < 0) {start = 0} #do not shift start.tmp if we are at the edge
  bin_start = round(start / bin.width + 1) #bin number of the first bin of orca matrix
  bin_end = bin_start + 249 #position of the last bin


  if (start / bin.width != start %/% bin.width) {
    message("Start/stop (ie ", start / 1e6, "Mb/", (start + scale)/1e6, "Mb) are not a multple of bin.width, rounded to ",
            (bin_start - 1) * bin.width/1e6, "Mb/", bin_end * bin.width / 1e6, "Mb.")
  }

  message("Creating matrix of ", nbins, "x", nbins, " bins with ", output, " counts of bins ", bin_start, " to ", bin_end,
          " at ", bin.width / 1e3, "kb resolution (i.e. from ", (bin_start - 1) * bin.width/1e6, "Mb to ", bin_end * bin.width / 1e6, "Mb).")

  #create empty matrix
  mat = Matrix::Matrix(0, nrow = nbins, ncol = nbins, sparse = TRUE)

  #read orca output
  if (output == "OE") {
    orca = utils::read.table(df_prediction.path, header = FALSE, sep = sep) %>% as.matrix() %>% exp
  }
  if (output == "Obs") {
    predictions = utils::read.table(df_prediction.path, header = FALSE, sep = sep) %>% as.matrix() %>% exp
    normats = utils::read.table(df_normmats.path, header = FALSE, sep = sep) %>% as.matrix()
    orca = predictions * normats
  }

  #add upper input in the full upper matrix
  mat[bin_start:bin_end, bin_start:bin_end][upper.tri(mat[bin_start:bin_end, bin_start:bin_end], diag = TRUE)] <- orca[upper.tri(orca, diag = TRUE)]

  return(methods::as(mat, "CsparseMatrix"))
}
