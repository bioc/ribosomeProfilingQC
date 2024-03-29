% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assignReadingFrame.R
\name{assignReadingFrame}
\alias{assignReadingFrame}
\title{Assign reading frame}
\usage{
assignReadingFrame(reads, CDS, txdb, ignore.seqlevelsStyle = FALSE)
}
\arguments{
\item{reads}{Output of \link{getPsiteCoordinates}}

\item{CDS}{Output of \link{prepareCDS}}

\item{txdb}{A TxDb object. If it is set, assign reading frame for all reads.
Default missing, only assign rading frame for reads in CDS.}

\item{ignore.seqlevelsStyle}{Ignore the sequence name style detection or not.}
}
\value{
An GRanges object of reads with reading frame information.
}
\description{
Set reading frame for each reads in CDS region to frame0,
frame1 and frame2.
}
\examples{
library(Rsamtools)
bamfilename <- system.file("extdata", "RPF.WT.1.bam",
                           package="ribosomeProfilingQC")
yieldSize <- 10000000
bamfile <- BamFile(bamfilename, yieldSize = yieldSize)
pc <- getPsiteCoordinates(bamfile, bestpsite=13)
pc.sub <- pc[pc$qwidth \%in\% c(29, 30)]
#library(GenomicFeatures)
library(BSgenome.Drerio.UCSC.danRer10)
#txdb <- makeTxDbFromGFF(system.file("extdata",
 #         "Danio_rerio.GRCz10.91.chr1.gtf.gz",
 #         package="ribosomeProfilingQC"),
 #         organism = "Danio rerio",
 #         chrominfo = seqinfo(Drerio)["chr1"],
 #         taxonomyId = 7955)
#CDS <- prepareCDS(txdb)
CDS <- readRDS(system.file("extdata", "CDS.rds",
                           package="ribosomeProfilingQC"))
pc.sub <- assignReadingFrame(pc.sub, CDS)
}
