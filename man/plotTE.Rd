% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTE.R
\name{plotTE}
\alias{plotTE}
\title{Plot translational efficiency}
\usage{
plotTE(
  TE,
  sample,
  xaxis = c("mRNA", "RPFs"),
  removeZero = TRUE,
  log2 = TRUE,
  theme = theme_classic(),
  type = "histogram",
  margins = "y",
  ...
)
}
\arguments{
\item{TE}{Output of \link{translationalEfficiency}}

\item{sample}{Sample names to plot.}

\item{xaxis}{What to plot for x-axis.}

\item{removeZero}{Remove the 0 values from plots.}

\item{log2}{Do log2 transform for TE or not.}

\item{theme}{Theme for ggplot2.}

\item{type, margins, ...}{Parameters pass to ggMarginal}
}
\value{
A ggExtraPlot object.
}
\description{
Scatterplot of RNA/RPFs level compared to the
translational efficiency.
}
\examples{
path <- system.file("extdata", package="ribosomeProfilingQC")
#RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#cnts <- countReads(RPFs, RNAs, gtf, level="gene")
cnts <- readRDS(file.path(path, "cnts.rds"))
fpkm <- getFPKM(cnts)
te <- translationalEfficiency(fpkm)
plotTE(te, 1)
}
