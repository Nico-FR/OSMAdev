#' @title Wild type metadata dataframe generation for specified regions
#'
#' @description
#' The main goal of this function is to create WT metatdata to run Orca 1Mb predictions only for the selected mutated windows.
#' The purpose of this function is to reduce the number of WT sequences needed to compare MT sequences. To do so, one WT sequence will be used as control for many MT sequences (number depending on how many mutation will be performed per mutated window).
#' For each 1Mb WT sequences, this function return the window that will be subject to mutations.
#' Those metadata will be used by createMetadaMT function to create the metadata of the mutant sequences (MT).
#'
#' @param DNAstring DNAString : the DNA sequence (e.g. chromosome).
#' @param DNAstring_name character : the sequence name (e.g 1 or chr1...)
#' @param start.windows numeric : vector of start positions of the windows that will contain the mutations.
#' @param stop.windows numeric : vector of stop positions of the windows that will contain the mutations.
#' @param max.window.size numeric : juxtaposed window will be reduce. If the width of the window is greater to max.window.size, the widow will be split. Must be smaller than 1Mb.
#' @param model_HFF <boolean> : use HFF training model for Orca predictions
#' @param model_ESC <boolean> : use ESC training model for Orca predictions
#'
#' @return dataframe
#'
#' @details
#' This function return a dataframe with metadatas relate to WildType sequences used to run Orca predictions :
#' Column description:
#' * chr <numeric> : the number of the chromosome
#' * start.window <numeric> : start of the window that will contain the mutations
#' * stop.window <numeric> : stop of the window that will contain the mutations
#' * start.mat <numeric> : the first position of the 1Mb sequence
#' * stop.mat <numeric> : the last position of the 1Mb sequence
#' * ID <character> : the name of the wild type sequence "WT_x.fa" with x an incremental number
#' * A,C,T,G <numeric> : 4 column with the frequences of each nucleotides on the 1Mb sequence.
#' * model = 1000000 : ORCA parameter cf ORCA documentation
#' * scale = 1000000 : ORCA parameter cf ORCA documentation
#' * wpos  500000 : ORCA parameter cf ORCA documentation
#' * mpos  500000 : ORCA parameter cf ORCA documentation
#' * model_HFF <0 or 1> : ORCA parameter. predict HFF matrices if 1
#' * model_ESC <0 or 1> : ORCA parameter. predict ESC matrices if 1
#'
#' @importFrom magrittr %>%
#' @importFrom Biostrings letterFrequency
#' @importFrom tibble tibble
#' @importFrom GenomicRanges GRanges start end tile reduce
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics unlist
#' @export


metadataWT_customWindow = function(DNAstring, DNAstring_name, start.windows, stop.windows, max.window.size, model_HFF, model_ESC){

  ####################################
  #testing parameters
  #DNAstring = BSgenome.Btaurus.UCSC.bosTau9::BSgenome.Btaurus.UCSC.bosTau9$chr1[1:2300000];DNAstring_name = "chr1";max.window.size = 256e3;model_HFF = TRUE;model_ESC = FALSE;gr1 = GenomicRanges::GRanges(seqnames = area_step2$chr, ranges = IRanges(area_step2$start.shuff, area_step2$stop.shuff), strand = "*");  start.windows = GenomicRanges::start(gr1);  stop.windows = GenomicRanges::end(gr1)
  ####################################

  if (max.window.size > 1000000){
    stop("window.size size must be smaller than 1Mb")
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
    start.mat = ifelse(start.window - 0.5e6 < 1, 1, start.window - 0.5e6),
    stop.mat = start.mat + 1e6 - 1,
    ID = paste0("WT_", 1:length(start.mat)),
    model = 1e6,
    scale = 1e6,
    wpos = 0.5e6,
    mpos = 0.5e6,
    model_HFF = ifelse(model_HFF, 1, 0),
    model_ESC = ifelse(model_ESC, 1, 0)
  )

  # Ensure that stop.mat does not exceed the length of the DNAstring
  metadataWT$stop.mat = pmin(metadataWT$stop.mat, dna.length)

  # Ensure that width is still 1Mb
  metadataWT$start.mat = metadataWT$stop.mat - 1e6 + 1

  ## computation of nucleotides frequences
  freq.WT = sapply(1:nrow(metadataWT),function(i){
    str = DNAstring[metadataWT$start.mat[i]:metadataWT$stop.mat[i]]
    Biostrings::letterFrequency(str, letters=c("A","T","G", "C"), as.prob = TRUE)
  }) %>% t

  ## add nucleotides frequencies to the metadata tibble
  metadataWT = cbind(metadataWT, freq.WT)

  # Return the dataframe
  return(metadataWT)
}
