#' @title DNA insertion
#'
#' @description insert of DNA segment into another DNA segment
#'
#' @param dna.string1 DNA input as DNAstring object to be inserted.
#' @param dna.string2 DNA input as DNAstring object (e.g chromosome sequence) where the insertion will take place.
#' @param ins.start position in bp of the insertion on DNA.string2
#' @param rev whether or not to reverse the inserted sequence, default = FALSE.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString
#' @export
#'
#' @examples
#' DNA.seq1 = Biostrings::DNAString("GATG")
#' DNA.seq2 = Biostrings::DNAString("AAAAACCCCCTTTTTAAAAA")
#'
#' #insertion startint at bp number 11:
#' insDNA2(dna.string1 = DNA.seq1, dna.string2 = DNA.seq2, ins.start = 11, rev = FALSE)
#'
#' #idem but reversed
#' insDNA2(dna.string1 = DNA.seq1, dna.string2 = DNA.seq2, ins.start = 11, rev = TRUE)
#'
insDNA2 <- function(dna.string1, dna.string2, ins.start, rev = FALSE) {

  mutated.seq = c(
    dna.string2[1:(ins.start - 1)],

    if (isTRUE(rev)) {rev(dna.string1)} else {dna.string1},

    dna.string2[(ins.start):length(dna.string2)])

  return(mutated.seq)
}

