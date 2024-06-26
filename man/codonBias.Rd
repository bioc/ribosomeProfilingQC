% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codonBias.R
\name{codonBias}
\alias{codonBias}
\title{Codon usage bias}
\usage{
codonBias(
  RPFs,
  gtf,
  genome,
  bestpsite = 13,
  readsLen = c(28, 29),
  anchor = "5end",
  ignore.seqlevelsStyle = FALSE,
  summary = TRUE,
  removeDuplicates = TRUE,
  ...
)
}
\arguments{
\item{RPFs}{Bam file names of RPFs.}

\item{gtf}{GTF file name for annotation or a TxDb object.}

\item{genome}{A BSgenome object.}

\item{bestpsite}{P site postion.}

\item{readsLen}{Reads length to keep.}

\item{anchor}{5end or 3end. Default is 5end.}

\item{ignore.seqlevelsStyle}{Ignore the sequence name style detection or not.}

\item{summary}{Return the summary of codon usage bias or full list.}

\item{removeDuplicates}{Remove the PCR duplicates or not. Default TRUE.}

\item{...}{Parameters pass to
\link[txdbmaker:makeTxDbFromGFF]{makeTxDbFromGFF}}
}
\value{
A list of data frame of codon count table if summary is TRUE.
 list 'reads' means the counts by raw reads.
 list 'reference' means the counts by sequence extracted from reference by
 the coordinates of mapped reads.
 Otherwise, return the counts (reads/reference) table for each reads.
}
\description{
Calculate the codon usage for the reads in the identified CDSs.
And then compared to the reference codon usage.
}
\examples{
path <- system.file("extdata", package="ribosomeProfilingQC")
RPFs <- dir(path, "RPF.*?\\\\.[12].bam$", full.names=TRUE)
gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
library(BSgenome.Drerio.UCSC.danRer10)
cb <- codonBias(RPFs[c(1,2)], gtf=gtf, genome=Drerio)
}
