#' @title DNA swapping within sequence
#'
#' @description From DNA sequence (e.g. chromosome), the function copy a DNA segment  which is moved to another position.
#' The positioning is done by overwriting the destination sequence of the same size.
#'
#' @details
#' Note that fo simplify, DNA positions are 0 base!
#' For example, if origin.start = 1e6 (1Mb) and width = 1e3, the first nucleotide position will be 1,000,001bp and the last 1,001,000.
#'
#' @param DNA.string DNA input as DNAstring object
#' @param center center of output DNA window in bp.
#' @param window output size in bp (centered around `center`). E.g, for +/-0.5 Mb on either side of `center`, use 1 Mb.
#' @param origin.start start of the segment to be swapped in bp.
#' @param width width of the segment to be swapped in bp.
#' @param target.start start of the segment to be overwriting in bp.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq = Biostrings::DNAString("AAAAAAAAAATTTTTTTTTTCCCCCCCCCC")
#' DNA.seq
#'
#' #swapping of the first 10 nucleotides (1:10) in the last 10 nucleotides (21:30)
#' #output window of 20 nucletodides (i.e. 11:30)
#' swapDNA(DNA.string = DNA.seq, center = 20, window = 20,
#'     origin.start = 0, width = 10, target.start = 20)
#'
swapDNA <- function(DNA.string, center, window, origin.start, width, target.start) {

  #window position
  w.start <- center - window %/% 2 + 1
  w.stop <- center + window %/% 2

  if (w.start < 0 | w.stop > length(DNA.string)) {
    stop("the window is outside the chr!")}

  if (w.start > target.start | w.stop < (target.start + width)) {
    stop("the target destination (start or stop) is outside the window!")}

  #sequence to swap
  seq.ori <- DNA.string[(origin.start + 1):(origin.start + width)]

  # window sequence
  seq <- DNA.string[w.start:w.stop]

  #sequence position to rm (relative to window)
  mut.start <- target.start - (center - window /2) + 1
  mut.stop <- target.start - (center - window /2) + width

  #swapping
  seq[mut.start:mut.stop] <- seq.ori

  return(seq)
}

