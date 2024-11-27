#' @title insert a DNA segment into another DNAsegment
#'
#' @description From two DNA sequence (e.g. chromosome), the function export a segment from the first DNA sequence (DNA.string1) into the second DNA sequence (DNA.string2).
#'
#' @details
#' By default the breakpoint of the output sequence is the middle of window. For example, if window is 10 bp, the breakpoint is between base pair number 5 and 6.
#'
#' @param DNA.string1 DNA input as DNAstring object to be inserted.
#' @param DNA.string2 DNA input as DNAstring object (e.g chromosome sequence) where the insertion will take place.
#' @param window output size in bp. E.g, for +/-0.5 Mb of each sequences, use 1 Mb.
#' @param breakpoint position in bp of the breakpoint (i.e position where the insertion occurs) on DNA.string2 sequence. If breakpoint = 11, the insertion is realised between bp number 10 and 11.
#' @param fpos final position in bp of the insertion according to the output sequence. By default fpos = "center" to center the insertion on the putût sequence.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq1 = Biostrings::DNAString("GG")
#' DNA.seq2 = Biostrings::DNAString("AAAAACCCCCTTTTTAAAAA")
#'
#' #output of 10bp with the insertion in the middle of the output:
#' insDNA(DNA.string1 = DNA.seq1, DNA.string2 = DNA.seq2, window = 10, breakpoint = 11, fpos = "center")
#'
#' #output of 10bp with the insertion at 6bp of output:
#' insDNA(DNA.string1 = DNA.seq1, DNA.string2 = DNA.seq2, window = 10, breakpoint = 11, fpos = 6)
#'
insDNA <- function(DNA.string1, DNA.string2, window, breakpoint, fpos = "center") {

  breakpoint <- breakpoint - 1

  if (fpos == "center") {

    fw = window - length(DNA.string1) # bp of DNA.string2 to keep

    if (fw %/% 2 == fw / 2) {
      origin.start1 = breakpoint - fw / 2 + 1
      start.seq <- DNA.string2[origin.start1:breakpoint]
      origin.stop2 = breakpoint + fw / 2
      end.seq <- DNA.string2[(breakpoint + 1):origin.stop2]
    } else {
      origin.start1 = breakpoint - fw %/% 2
      start.seq <- DNA.string2[origin.start1:breakpoint]
      origin.stop2 = breakpoint + fw %/% 2
      end.seq <- DNA.string2[(breakpoint + 1):origin.stop2]
    }

    message(paste0("Insertion start (relative to output sequence) is bp number ", length(start.seq) + 1))
  }

  if (is.numeric(fpos)) {
    if (fpos > (breakpoint + 1)) {
      stop("fpos must be <= to breakpoint!")
    }
    if (fpos > window) {
      stop("fpos must be <= to window!")
    }

    origin.start1 = breakpoint - fpos + 2
    start.seq <- DNA.string2[origin.start1:breakpoint]
    if (origin.start1 > breakpoint) {start.seq <- Biostrings::DNAString("")}

    origin.stop2 = breakpoint + (window - length(start.seq)) - length(DNA.string1)
    if (origin.stop2 > length(DNA.string2)) {stop("length of DNA.string2 is insufficient!")}
    end.seq <- DNA.string2[(breakpoint + 1):origin.stop2]
    if (origin.stop2 < breakpoint) {end.seq <- Biostrings::DNAString("")}
  }

  #output sequence
  seq = c(start.seq, DNA.string1, end.seq)
  return(seq[1:window])
}
