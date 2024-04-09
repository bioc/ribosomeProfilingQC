#' Plot translational efficiency
#' @description Scatterplot of RNA/RPFs level compared to the
#' translational efficiency.
#' @param TE Output of \link{translationalEfficiency}
#' @param sample Sample names to plot.
#' @param xaxis What to plot for x-axis.
#' @param removeZero Remove the 0 values from plots.
#' @param log2 Do log2 transform or not.
#' @param type,margins,... Parameters pass to ggMarginal
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot geom_point geom_text .data facet_wrap
#' @importFrom ggExtra ggMarginal
#' @importFrom stats as.formula coef na.omit
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' #RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' #RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' #gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' #cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' fpkm <- getFPKM(cnts)
#' te <- translationalEfficiency(fpkm)
#' plotTE(te, 1)

plotTE <- function(TE, sample, xaxis=c("mRNA", "RPFs"),
                   removeZero=TRUE, log2=TRUE,
                   type = 'histogram', margins = 'y',
                   ...){
  if(!is.list(TE)){
    stop("TE must be output of translationalEfficiency.")
  }
  if(!any(c("RPFs", "mRNA", "TE") %in% names(TE))){
    stop("TE must be output of translationalEfficiency.")
  }
  xaxis <- match.arg(xaxis)
  mRNA <- TE$mRNA
  RPFs <- TE$RPFs
  x <- TE[[xaxis]]
  TE <- TE$TE
  if(missing(sample)){
    sample <- colnames(TE)
  }
  if(!is.numeric(sample)){
    sample <- which(colnames(TE) %in% sample)
  }
  single_sample <- length(sample)==1
  df <- lapply(sample, function(.sample){
    data.frame(x = x[, .sample],
               TE = TE[, .sample],
               sample = .sample,
               RPFs = RPFs[, .sample],
               mRNA = mRNA[, .sample])
  })
  
  df <- do.call(rbind, df)
  if(removeZero){
    keep <- df$RPFs>0 & df$mRNA>0
    if(sum(keep)<1){
      stop("No data available for plotting.")
    }
    df <- df[keep, , drop=FALSE]
  }
  df$sample <- colnames(TE)[df$sample]
  xlab <- paste(xaxis, "level")
  ylab <- "Translational Efficiency"
  if(log2){
    df$x <- log2(df$x)
    df$TE <- log2(df$TE)
    xlab <- paste('log2 transformed', xlab)
    ylab <- paste('log2 transformed', ylab)
  }

  lm_eqn <- function(df){
    df <- df[!is.infinite(df$x) & !is.infinite(df$TE), , drop=FALSE]
    df <- split(df, df$sample)
    out <- lapply(df, function(.df){
      m <- lm(as.formula("TE ~ x"), .df, na.action = na.omit)
      eq <- substitute(italic("TE") == a + b %.% italic(x)*","~~italic("r")^2~"="~r2, 
                       list(a = format(unname(coef(m)[1]), digits = 2),
                            b = format(unname(coef(m)[2]), digits = 2),
                            r2 = format(summary(m)$r.squared, digits = 3),
                            x = xaxis))
      as.character(as.expression(eq))
    })
    out <- data.frame(sample=names(out), label=unlist(out))
  }
  p <- ggplot(df, aes(.data$x, .data$TE)) + geom_point() + 
    xlab(xlab) + ylab(ylab) +
    geom_smooth(method="lm", se=TRUE)
  p <- p + geom_text(data=lm_eqn(df),
                     aes(x=-Inf, y = Inf, label=.data$label),
                     hjust = 0, vjust = 1,
                     parse = TRUE)
  if(single_sample){
    p <- ggMarginal(p, type=type, margins=margins, ...)
  }else{
    p <- p +
      facet_wrap(facets = as.formula('~sample'))
  }
  p
}
