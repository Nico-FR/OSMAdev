#' Wrapper for Orca data preparation
#'
#' @description description
#' This function prepares the data for Orca by generating WildType and Mutated metadata, writing the corresponding fasta files, and saving the metadata as tsv files.
#' This function run the following functions:
#' - metadataWT_customWindow or metadataWT_seqWide to generate WildType metadatas
#' - metadataMT to generate Mutated metadatas
#' - writeFastaWT to write WildType fasta files
#' - writeFastaMT to write Mutated fasta files
#'
#'
#' @inheritParams metadataWT_customWindow
#' @inheritParams metadataWT_seqWide
#' @inheritParams metadataMT
#' @inheritParams writeFastaMT
#' @inheritParams writeFastaWT
#'
#' @returns <message> : a message indicating the number of fasta files saved and the output directory.
#'
#' @export
#'
#'
writeOrcaData <- function(DNAstring, DNAstring_name, window.size = NULL, start.windows = NULL, stop.windows = NULL, max.window.size = NULL, model_HFF, model_ESC, mutated.width, workdir, gzip = FALSE, rep = 1) {

  if (!is.null(max.window.size) && (!is.null(start.windows) && !is.null(stop.windows) && is.null(window.size))) { # if max.window.size is provided with start and stop windows (i.e run metadataWT_customWindow function)

    #check if windows are multiple of mutated.width
    if ((stop.windows - start.windows + 1) %% mutated.width %>% mean() != 0) {
      stop("window sizes (i.e: start.windows - stop.windows + 1) must be multiple of mutated.width")
    }

    # create metadata for custom windows
    metadataWT = metadataWT_customWindow(DNAstring = DNAstring, DNAstring_name = DNAstring_name, start.windows = start.windows,
                                         stop.windows = stop.windows, max.window.size = max.window.size, model_HFF = model_HFF, model_ESC = model_ESC)

  } else if (is.null(max.window.size) && (is.null(start.windows) && is.null(stop.windows) && !is.null(window.size))) { # if only window.size is provided (i.e run metadataWT_seqWide function)

    #check if window.size are multiple of mutated.width
    if (window.size %% mutated.width != 0) {
      stop("window.size must be multiple of mutated.width")
    }

    # create metadata for wide sequence
    metadataWT = metadataWT_seqWide(DNAstring = DNAstring, DNAstring_name = DNAstring_name, window.size = window.size, model_HFF = model_HFF, model_ESC = model_ESC)

  } else {

    # otherwise stop
    stop("max.window.size, start.windows and stop.windows must be provided together. Otherwise only window.size must be provided.")
  }

  metadataMT = metadataMT(metadataWT, mutated.width, rep = rep)

  #write the fasta files
  writeFastaWT(DNAstring = DNAstring, metadataWT = metadataWT, workdir = workdir, gzip = gzip)
  writeFastaMT(DNAstring = DNAstring, metadataWT = metadataWT, metadataMT = metadataMT, workdir = workdir, gzip = gzip)

  #write the metadata files
  options(scipen = 100)
  utils::write.table(metadataWT, file=paste0(workdir, "metadataWT.tsv"), quote=FALSE, sep='\t', row.names = F)
  utils::write.table(metadataMT, file=paste0(workdir, "metadataMT.tsv"), quote=FALSE, sep='\t', row.names = F)
  options(scipen = 0)

  message("Metadatas saved in ", workdir)
}











