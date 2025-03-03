% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insDNA.R
\name{insDNA}
\alias{insDNA}
\title{DNA insertion}
\usage{
insDNA(dna.string1, dna.string2, start, rev = FALSE)
}
\arguments{
\item{dna.string1}{DNA input as DNAstring object to be inserted.}

\item{dna.string2}{DNA input as DNAstring object (e.g chromosome sequence) where the insertion will take place.}

\item{start}{position in bp of the insertion on DNA.string2}

\item{rev}{whether or not to reverse the inserted sequence, default = FALSE.}
}
\value{
DNAstring
}
\description{
insert of DNA segment into another DNA segment
}
\examples{
DNA.seq1 = Biostrings::DNAString("GATG")
DNA.seq2 = Biostrings::DNAString("AAAAACCCCCTTTTTAAAAA")

# insertion starting at bp number 11:
insDNA(dna.string1 = DNA.seq1, dna.string2 = DNA.seq2, start = 11, rev = FALSE)

# idem but reversed
insDNA(dna.string1 = DNA.seq1, dna.string2 = DNA.seq2, start = 11, rev = TRUE)

}
