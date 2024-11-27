#' @title create reciprocal translocation between two sequences
#'
#' @description From two DNA sequence (e.g. chromosome), the function fused a segment of the 2 sequences to create a new one.
#'
#' @details
#' By default the breakpoint of the output sequence is the middle of window. For example, if window is 10 bp, the breakpoint is between base pair number 5 and 6.
#'
#' @param DNA.string1 DNA input as DNAstring object.
#' @param DNA.string2 DNA input as DNAstring object.
#' @param window output size in bp. E.g, for +/-0.5 Mb of each sequences, use 1 Mb.
#' @param breakpoint position in bp of the breakpoint (i.e position where the change occurs) according to the exported sequence. Default is NULL to use the center of the window as breakpoint (see examples).
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq1 = Biostrings::DNAString("AAAAAGGGGG")
#' DNA.seq2 = Biostrings::DNAString("CCCCCTTTTT")
#'
#' # reciprocal translocation of the DNA.seq1 with DNA.seq2
#' trcpDNA(DNA.string1 = DNA.seq1, DNA.string2 = DNA.seq2, window = 10)
#'
#' # reciprocal translocation of the DNA.seq1 with DNA.seq2,
#' # which sequence change at third bp (breakpoint = 3)
#' trcpDNA(DNA.string1 = DNA.seq1, DNA.string2 = DNA.seq2, window = 10, breakpoint = 3)
#'
trcpDNA <- function(DNA.string1, DNA.string2, window, breakpoint = NULL) {

  if (is.null(breakpoint)) {breakpoint = window %/% 2 + 1}
  if (window < breakpoint) {stop("the breakpoint position is outside the window!")}

  #sequences to be merged
  origin.start1 <- length(DNA.string1) - breakpoint + 2
  origin.stop2 <- window - breakpoint + 1

  #sanity check
  if (origin.start1 >= length(DNA.string1)) {stop("the length of DNA.string1 is insufficient!")}
  if (origin.stop2 >= length(DNA.string2)) {stop("the length of DNA.string2 is insufficient!")}

  #begin of sequence
  seq = c(DNA.string1[origin.start1:length(DNA.string1)], DNA.string2[1:origin.stop2])

  return(seq)
}
