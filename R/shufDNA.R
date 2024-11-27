#' @title DNA shuffle
#'
#' @description From DNA sequence (e.g. chromosome), the function generate "a random" sequence of a DNA segment. The random sequence is generated according to the nucleotide frequencies of the supplied DNA set.
#'
#' @details
#' Note that for simplify, DNA positions are 0 base!
#' For example, if shuffle.start = 1e6 (1Mb) and width = 1e3, the first nucleotide position will be 1,000,001bp and the last 1,001,000.
#'
#' @param DNA.string DNA input as DNAstring object
#' @param center center of output DNA window in bp.
#' @param window output size in bp (centered around `center`). For +/-0.5 Mb on either side of `center`, use 1 Mb.
#' @param shuffle.start start of the segment to shuffle in bp.
#' @param width width of the segment to shuffle in bp.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString alphabetFrequency
#' @export
#'
#' @examples
#' DNA.seq = Biostrings::DNAString("AAAAAAAAAATTTTTTTTTTCCCCCCCCCC")
#' DNA.seq
#'
#' #shuffle of the last 10 nucleotides (21:30)
#' #output window of 30 nucletodides (i.e. 1:30)
#' shufDNA(DNA.string = DNA.seq, center = 15, window = 30,
#'     shuffle.start = 20, width = 10)
#'
shufDNA <- function(DNA.string, center, window, shuffle.start, width) {

  #window position
  w.start <- center - window %/% 2 + 1
  w.stop <- center + window %/% 2

  if (w.start < 0) {w.start = 1 ; w.stop = window}
  if (w.stop > length(DNA.string)) {
    stop("the window is outside the DNA set!")}

  if (w.start > shuffle.start | w.stop < (shuffle.start + width)) {
    stop("the segment to shuffle (start or stop) is outside the window!")}

  #freq of window
  seq = DNA.string[w.start:w.stop]
  freqs = Biostrings::alphabetFrequency(seq, baseOnly=TRUE, as.prob = TRUE)

  mut.start = shuffle.start - (center - window / 2) + 1
  mut.stop = shuffle.start - (center - window / 2) + width

  # shuffle
  mut <- sample(names(freqs)[names(freqs) %in% c("A","T","C","G")],
                width, replace = TRUE,
                prob = freqs[names(freqs) %in% c("A","T","C","G")])
  mut <- DNAString(paste(mut, collapse = ""))

  #overwrite muated sequence
  seq[mut.start:mut.stop] <- mut

  return(seq)
}

