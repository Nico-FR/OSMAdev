#' @title reverse a DNA segment (inversion)
#'
#' @description to do
#'
#'
#' @param dna.string DNA input as DNAstring object (e.g chromosome sequence)
#' @param inv.start position in bp of the first nucleotide to reverse
#' @param inv.stopt position in bp of the last nucleotide to reverse
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq = Biostrings::DNAString("AAAAACCCCCTTTTTAAAAA")
#'
#' #inversion of CCCCCTTTTT nucleotides:
#' invDNA(dna.string = DNA.seq, inv.start = 6, inv.stop = 15)
#'
invDNA <- function(dna.string, inv.start, inv.stop) {

  mutated.seq = c(

    dna.string[1:(inv.start - 1)],

    rev(dna.string[inv.start:inv.stop]),

    dna.string[(inv.stop + 1):length(dna.string)])

  return(mutated.seq)
}
