#' @title Wild type metadata dataframe generation sequence wide.
#'
#' @description
#' The main goal of this function is to create WT metatdata to run Orca 1Mb predictions all along the sequence provided.
#' The purpose of this function is to reduce the number of WT sequences needed to compare MT sequences. To do so, one WT sequence will be used as control for many MT sequences (number depending on how many mutation will be performed per window.size).
#' For each 1Mb WT sequences, this function return the window that will be subject to mutations, ensuring that these mutated windows are adjacent along the sequence
#' Those metadata will be used by createMetadaMT function to create the metadata of the mutant sequences (MT).
#'
#' @param DNAstring DNAString : the DNA sequence (e.g. chromosome).
#' @param DNAstring_name character : the sequence name (e.g "chr1")
#' @param window.size numeric : size of adjacent windows in bp which indicate mutations area at the center of the 1Mb sequences. Must be smaller than 1Mb (see details)
#' @param model_HFF <boolean> : use HFF training model for Orca predictions
#' @param model_ESC <boolean> : use ESC training model for Orca predictions
#'
#' @return dataframe
#'
#' @details
#' The window.size parameter is used to define the start and stop of the mutated area for each 1Mb sequences. This window is positioned at the center of the 1Mb sequences that will be predicted.
#' For example: For the first 1Mb sequence (1 to 1e6 bp), if window = 256000, the fist mutation will start at 372001 (1e6 / 2 - 256e3 / 2 + 1) and the last mutation will start at 628000 (1e6 / 2 - 256e3 / 2).
#'
#' This function return a dataframe with metadatas relate to WildType sequences used to run Orca predictions :
#' Column description:
#' * chr <numeric> : the name of the sequence
#' * start.window <numeric> : start of the window that will contain the mutations
#' * stop.window <numeric> : stop of the window that will contain the mutations
#' * start.mat <numeric> : the first position of the 1Mb sequence
#' * stop.mat <numeric> : the last position of the 1Mb sequence
#' * ID <character> : the name of the wild type sequence "WT_x.fa" with x an incremental number
#' * A,C,T,G <numeric> : 4 column with the frequencies of each polynucleotides on the 1Mb sequence.
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
#' @export


metadataWT_seqWide = function(DNAstring, DNAstring_name, window.size, model_HFF, model_ESC){

  ####################################
  #testing parameters
  #DNAstring = BSgenome.Btaurus.UCSC.bosTau9::BSgenome.Btaurus.UCSC.bosTau9$chr1[1:2534110];DNAstring_name = "chr1";window.size = 256e3;model_HFF = TRUE;model_ESC = FALSE
  ####################################

  if (window.size > 1000000){
    stop("window.size size must be smaller than 1Mb")
  }

  dna.length = length(DNAstring)
  first.start.window = 0.5e6 - window.size / 2 + 1 #start of the first mutation
  last.start.window = length(DNAstring) - first.start.window + 1 - window.size #stop of the last mutation


  ### create metadata
  metadataWT = tibble::tibble(
    chr = DNAstring_name,
    start.window = seq(first.start.window, length(DNAstring) - first.start.window + 1 - window.size, by = window.size),
    stop.window = start.window + window.size -1,
    start.mat = start.window - first.start.window + 1,
    stop.mat = start.mat + 1e6 - 1,
    ID = paste0("WT_", 1:length(start.window)),
    model = 1e6,
    scale = 1e6,
    wpos = 0.5e6,
    mpos = 0.5e6,
    model_HFF = ifelse(model_HFF, 1, 0),
    model_ESC = ifelse(model_ESC, 1, 0)
  )

  ## compuation of nucleotides frequences
  freq.WT = sapply(1:nrow(metadataWT),function(i){
    str = DNAstring[metadataWT$start.mat[i]:metadataWT$stop.mat[i]]
    Biostrings::letterFrequency(str, letters=c("A","T","G", "C"), as.prob = TRUE)
  }) %>% t

  ## add nucleotides frequencies to the metadata tibble
  metadataWT = cbind(metadataWT, freq.WT)

  # Return the dataframe
  return(metadataWT)
}
