#' Normalize the TE by Loess
#' @description
#' Fitting the translational efficiency values with the mRNA value by \link[stats]{loess}.
#' @param TE output of \link{translationalEfficiency}.
#' @param log2 logical(1L). Do log2 transform for TE or not. If TE value is not
#' log2 transformed, please set it as TRUE.
#' @param pseudocount The number will be add to sum of reads count to avoid X/0.
#' @param span,family.loess Parameters will be passed to \link[stats]{loess}
#' @return A list with RPFs, mRNA levels and TE as a matrix with
#' log2 transformed translational efficiency.
#' @importFrom stats loess predict
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' fpkm <- getFPKM(cnts)
#' te <- translationalEfficiency(fpkm)
#' te1 <- normalizeTEbyLoess(te)
#' plotTE(te)
#' plotTE(te1, log2=FALSE)
normalizeTEbyLoess <- function(TE, log2=TRUE,
                               pseudocount=1e-3, span=2/3,
                               family.loess="symmetric"){
  if(!is.list(TE)){
    stop("TE must be output of translationalEfficiency.")
  }
  if(!any(c("RPFs", "mRNA", "TE") %in% names(TE))){
    stop("TE must be output of translationalEfficiency.")
  }
  stopifnot(is.logical(log2))
  log2 <- log2[1]
  for(j in seq.int(ncol(TE[["TE"]]))){
    x <- TE[['mRNA']][, j, drop=TRUE]
    y <- TE[["TE"]][, j, drop=TRUE]
    x <- log2(x+pseudocount)
    if(log2){
      y <- log2(y+pseudocount)
    }
    index <- order(x)
    xx <- x[index]
    yy <- y[index]
    aux <- loess(as.formula('yy ~ xx'), span=span,
                 degree=1, family=family.loess)
    aux <- predict(aux, data.frame(xx=x))
    TE[["TE"]][, j] <- y-aux
  }
  TE
}