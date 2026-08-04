#' @title Wild type metadata dataframe generation for specified regions with AlphaGenome
#'
#' @description
#' The main goal of this function is to create WT metatdata to run AlphaGenome 1048576 bp predictions only for the selected mutated windows.
#' Those metadata will be used by createMetadaMT function to create the metadata of the mutant sequences (MT).
#'
#' @param DNAstring DNAString : the DNA sequence (e.g. chromosome).
#' @param DNAstring_name character : the sequence name (e.g 1 or chr1...)
#' @param start.windows numeric : vector of start positions of the windows that will contain the mutations.
#' @param stop.windows numeric : vector of stop positions of the windows that will contain the mutations.
#' @param max.window.size numeric : juxtaposed window will be reduce. If the width of the window is greater to max.window.size, the widow will be split. Must be smaller than 1048576.
#' @param organism character : organism for AlphaGenome prediction (e.g., "human" or "mouse")
#' @param ontology_term character : ontology term for AlphaGenome prediction
#'
#' @return dataframe
#'
#' @details
#' This function return a dataframe with metadatas relate to WildType sequences used to run AlphaGenome predictions :
#' Column description:
#' * chr <numeric> : the number of the chromosome
#' * start.window <numeric> : start of the window that will contain the mutations
#' * stop.window <numeric> : stop of the window that will contain the mutations
#' * start.pred <numeric> : the first position of the 1048576 bp sequence
#' * stop.pred <numeric> : the last position of the 1048576 bp sequence
#' * ID <character> : the name of the wild type sequence "WT_x.fa" with x an incremental number
#' * A,C,T,G <numeric> : 4 column with the frequences of each nucleotides on the 1048576 bp sequence.
#' * scale = 1048576 : ALPHAGENOME parameter
#' * organism <character> : ALPHAGENOME parameter
#' * ontology_term <character> : ALPHAGENOME parameter
#'
#' @importFrom magrittr %>%
#' @importFrom Biostrings letterFrequency
#' @importFrom tibble tibble
#' @importFrom GenomicRanges GRanges start end tile reduce
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics unlist
#' @export


metadataWT_customWindow_Alphagenome = function(DNAstring, DNAstring_name, start.windows, stop.windows, max.window.size, organism, ontology_term){

  if (max.window.size > 1048576){
    stop("window.size size must be smaller than 1048576")
  }

  if (length(start.windows) != length(stop.windows)) {
    stop("start.windows and stop.windows must have the same length")
  }

  dna.length = length(DNAstring)

  #create GRange and split ranges > max.window.size
  gr1 = GenomicRanges::GRanges(
    seqnames = DNAstring_name,
    ranges = IRanges::IRanges(start = start.windows, end = stop.windows),
    strand = "*"
  ) %>% GenomicRanges::reduce() %>% GenomicRanges::tile(width = max.window.size) %>% BiocGenerics::unlist()

  ### create metadata
  metadataWT = tibble::tibble(
    chr = DNAstring_name,
    start.window = GenomicRanges::start(gr1),
    stop.window = GenomicRanges::end(gr1),
    start.pred = ifelse(start.window - 524288 < 1, 1, start.window - 524288),
    stop.pred = start.pred + 1048576 - 1,
    ID = paste0("WT_", 1:length(start.pred)),
    scale = 1048576,
    organism = organism,
    ontology_term = ontology_term
  )

  # Ensure that stop.pred does not exceed the length of the DNAstring
  metadataWT$stop.pred = pmin(metadataWT$stop.pred, dna.length)

  # Ensure that width is still 1048576
  metadataWT$start.pred = metadataWT$stop.pred - 1048576 + 1

  ## computation of nucleotides frequences
  freq.WT = sapply(1:nrow(metadataWT),function(i){
    str = DNAstring[metadataWT$start.pred[i]:metadataWT$stop.pred[i]]
    Biostrings::letterFrequency(str, letters=c("A","T","G", "C"), as.prob = TRUE)
  }) %>% t

  ## add nucleotides frequencies to the metadata tibble
  metadataWT = cbind(metadataWT, freq.WT)

  # Return the dataframe
  return(metadataWT)
}
