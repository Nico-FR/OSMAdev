#' @title create copy number variations of sequences
#'
#' @description From DNA sequence (e.g. chromosome), the function create X copy of a DNA segment (e.g. deletion, duplication...).
#'
#'
#' @param dna.string DNA input as DNAstring object (e.g chromosome sequence)
#' @param cnv.stop stop of the CNV in bp.
#' @param cnv.start start of the CNV in bp.
#' @param nb_copy final number of copy of the CNV (0 copy for deletion).
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
#' cnvDNA2(dna.string = DNA.seq, cnv.start = 16, cnv.stop = 17, nb_copy = 2)
#'
#' #deletion of "T" nucleotides (16:17)
#' cnvDNA2(dna.string = DNA.seq, cnv.start = 16, cnv.stop = 17, nb_copy = 0)
#'
cnvDNA2 <- function(dna.string, cnv.start, cnv.stop, nb_copy) {

  mutated.seq = c(
    dna.string[1:(cnv.start - 1)],
    rep(dna.string[cnv.start:cnv.stop], nb_copy),
    dna.string[(cnv.stop + 1):length(dna.string)])

  return(mutated.seq)
}
