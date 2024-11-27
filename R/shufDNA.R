#' @title DNA shuffle
#'
#' @description Shuffle of a DNA segment with random or custom nucleotides calling.
#'
#' @details The nucleotides calling can be based with differents nucleotides probabilities (see examples):
#'  1) the nucleotide frequency of the full sequence provided,
#'  2) the nucleotide frequency of the segment to shuffle,
#'  3) an equivalent frequency for each nucleotide (i.e random calling).
#'  4) custom probabilities of each nucleotides.
#'
#' @param dna.string DNA input as DNAstring object
#' @param start start of the segment to shuffle in bp.
#' @param stop stop of the segment to shuffle in bp.
#' @param probability A,T,G and C nucleotides probabilities respectively. Default is c(0.25, 0.25, 0.25, 0.25) to randomly call each nucleotide. Use probability = "full" to perform randomization based on the nucleotide probability of the entire sequence supplied. Use the “partial” parameter to take into account only the probability of the segment to be randomized.
#'
#' @return DNAstring
#'
#' @importFrom Biostrings DNAString letterFrequency
#' @export
#'
#' @examples
#' dna.seq = Biostrings::DNAString("AAAAATTTTTGGGGGAAAAA")
#' dna.seq
#'
#' # shuffle of T and G nucleotides (6:15) with random calling of nucleotides,
#' # i.e probability = (0.25, 0.25, 0.25, 0.25):
#' shufDNA(dna.string = dna.seq, start = 6, stop = 15)
#'
#' # shuffle of T and G nucleotides (6:15) based on the probability of the entire sequence provided,
#' # i.e equivalent to probability = c(0.5, 0.25, 0.25, 0) in this example::
#' shufDNA(dna.string = dna.seq, start = 6, stop = 15, probability = "full")
#'
#' # idem but based on the probability of the segment to shuffle,
#' # i.e equivalent to probability = c(0, 0.5, 0.5, 0) in this example:
#' shufDNA(dna.string = dna.seq, start = 6, stop = 15, probability = "partial")
#'
#' # idem but with custom probabilities:
#' # probability = c(0.5, 0, 0, 0.5) for A,T,G and C respectively:
#' shufDNA(dna.string = dna.seq, start = 6, stop = 15, probability = c(0.5, 0, 0, 0.5))
#'
shufDNA <- function(dna.string, start, stop, probability = c(0.25, 0.25, 0.25, 0.25)) {

  shuffle.width = stop - start + 1

  if (unique(probability == "full")) {
    probability = Biostrings::letterFrequency(dna.string, letters=c("A","T","G", "C"), as.prob = TRUE)
    }

  if (unique(probability == "partial")) {
    probability = Biostrings::letterFrequency(dna.string[start:stop], letters=c("A","T","G", "C"), as.prob = TRUE)
  }

  #calling/shuffle
  mut <- sample(c("A","T","G","C"),
                shuffle.width, replace = TRUE,
                prob = probability)

  mut <- DNAString(paste(mut, collapse = ""))

  #overwrite muated sequence
  dna.string[start:stop] <- mut

  return(dna.string)
}

