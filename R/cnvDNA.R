#' @title create copy number variations of sequences
#'
#' @description From DNA sequence (e.g. chromosome), the function create X copy of a DNA segment (e.g. deletion, duplication...).
#'
#'
#' @param DNA.string DNA input as DNAstring object (e.g chromosome sequence)
#' @param window output size in bp.
#' @param cnv.start start of the CNV in bp.
#' @param width width of the CNV in bp.
#' @param fpos final position in bp of the CNV according to the output sequence. By default fpos = "center" to center the CNV on the output sequence.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq = Biostrings::DNAString("AAAAAAAAAAAAAAGTTGCCCCCCCCCCCCCC")
#' DNA.seq
#'
#' #duplication of "T" nucleotides (16:17)
#' #output window of 20 nucletodides with de duplication start at the center (i.e 11th bp)
#' cnvDNA(DNA.string = DNA.seq, window = 20,
#'     cnv.start = 16, width = 2, nb_copy = 2, fpos = "center")
#'
#' #output window of 20 nucletodides with de duplication start at 6bp
#' cnvDNA(DNA.string = DNA.seq, window = 20,
#'     cnv.start = 16, width = 2, nb_copy = 2, fpos = 6)
#'
#' #deletion of "T" nucleotides (16:17)
#' cnvDNA(DNA.string = DNA.seq, window = 20,
#'     cnv.start = 16, width = 2, nb_copy = 0)
#'
cnvDNA <- function(DNA.string, window, cnv.start, width, nb_copy, fpos = "center") {

  if (fpos == "center") {
    fpos = window %/% 2 + 1 # final position of CNV
  }

  if (fpos > window) {
    stop("fpos must be <= to window!")
  }

  #before CNV
  origin.stop1 = cnv.start - 1
  origin.start1 = cnv.start - fpos + 1
  start.seq <- DNA.string[origin.start1:origin.stop1]
  if (origin.stop1 < origin.start1) {
    origin.start1 = origin.stop1 <- NULL
    start.seq <- Biostrings::DNAString("")}

  #CNV
  cnv.seq = rep(DNA.string[cnv.start:(cnv.start + width - 1)], nb_copy)

  #after CNV
  origin.start2 = cnv.start + width
  origin.stop2 = origin.start2 + window - length(cnv.seq) - length(start.seq)
  if (origin.stop2 > length(DNA.string)) {stop("length of DNA.string is insufficient!")}
  if (origin.stop2 < origin.start2) {
    end.seq <- Biostrings::DNAString("")
    origin.stop2 = origin.start2 <- NULL
    } else {
      end.seq <- DNA.string[origin.start2:origin.stop2]}

  message(
    paste0("1st nucleotide position = ", min(c(origin.start1, origin.stop1, origin.start2, origin.stop2, cnv.start)))
    )
  seq = c(start.seq, cnv.seq, end.seq)
  return(seq[1:window])

}

