#' @title reverse a DNA segment (inversion)
#'
#' @description to do
#'
#'
#' @param dna.string DNA input as DNAstring object (e.g chromosome sequence)
#' @param start position in bp of the first nucleotide to reverse
#' @param stop position in bp of the last nucleotide to reverse
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString reverseComplement
#' @export
#'
#' @examples
#' DNA.seq = Biostrings::DNAString("AAAAACCCCCAAAAACCCCC")
#'
#' # inversion of CCCCCAAAAA nucleotides:
#' invDNA(dna.string = DNA.seq, start = 6, stop = 15)
#'
invDNA <- function(dna.string, start, stop) {

  mutated.seq = c(

    dna.string[1:(start - 1)],

    Biostrings::reverseComplement(dna.string[start:stop]),

    dna.string[(stop + 1):length(dna.string)])

  return(mutated.seq)
}
