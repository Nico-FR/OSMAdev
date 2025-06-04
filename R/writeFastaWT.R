#' @title metadata WT to fasta sequences files
#'
#' @param DNAstring <DNAString> : the sequence of the chromosome, it must be the same sequence that was use to generate the WildType metadata dataFrame (metadataWT_seqWide or metadataWT_customWindow).
#' @param metadataWT <dataFrame> : a dataframe generate by metadataWT_seqWide or metadataWT_customWindow.
#' @param workdir <character> : the working directory where the fasta files will be saved. Default is "./".
#' @param gzip <logical> : if TRUE, the fasta files will be saved in gzip format. Default is TRUE.
#'
#' @return  <message> : a message indicating the number of fasta files saved and the output directory.
#'
#' @importFrom Biostrings writeXStringSet DNAStringSet
#' @importFrom progress progress_bar
#' @export


writeFastaWT = function(DNAstring, metadataWT, workdir ="./", gzip = TRUE) {

  ##############################
  # testing parameters
  #DNAstring = BSgenome.Btaurus.UCSC.bosTau9::BSgenome.Btaurus.UCSC.bosTau9$chr1[1:2534110]; gz = TRUE; workdir="/media/nmary/DONNEES/bureau_nico/Recherche/R/NoteBooks/"
  ###############################

  #create the working directory if it doesn't exist
  if (!dir.exists(workdir)) {
    dir.create(workdir, recursive = TRUE)
  }

  #create the output if it doesn't exist
  outdir = paste0(workdir, "/fasta/")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  #create DNAStringSet
  WT.fa.lst = Biostrings::DNAStringSet(DNAstring, start = metadataWT$start.mat, end = metadataWT$stop.mat)
  names(WT.fa.lst) = metadataWT$ID
  message("Writing ", nrow(metadataWT), " WT fasta files...")

  #initiate progress bar
  pb <- progress::progress_bar$new(total = length(WT.fa.lst),
                                   format = "[:bar] :percent in :elapsed",
                                   clear = FALSE, width= 60)
  pb$tick(0)

  #save the fasta file
  for (i in 1:length(WT.fa.lst)) {

    if (isFALSE(gzip)) {
      #write the fasta file
      Biostrings::writeXStringSet(WT.fa.lst[i], filepath = paste0(outdir, names(WT.fa.lst[i]), ".fa"))
    }

    if (isTRUE(gzip)) {
      #write the fasta file
      Biostrings::writeXStringSet(WT.fa.lst[i], filepath = paste0(outdir, names(WT.fa.lst[i]), ".fa.gzip"), compress = TRUE)
    }
    pb$tick()
  }

  return(message(nrow(metadataWT), " fasta files saved in ", outdir))
}

