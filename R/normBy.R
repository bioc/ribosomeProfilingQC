#' Normalization by edgeR, DESeq2 or RUVSeq
#' @description Normalization by multiple known methods
#' @param counts Output of \link{countReads}
#' @param method Character(1L) to indicate the method for normalization.
#' @param ... parameters will be passed to \link{normByRUVs} or \link{getFPKM}
#' @return Normalized counts list
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' cnts <- readRDS(file.path(path, "cnts.rds"))
#' norm <- normBy(cnts, method = 'edgeR')
#' norm2 <- normBy(cnts, method = 'DESeq2')
#' 
normBy <- function(counts, method = c('edgeR', 'DESeq2', 'RUVs', 'fpkm'),
                   ...){
  if(!any(c("RPFs", "mRNA") %in% names(counts))){
    stop("counts must be output of coutReads.")
  }
  method <- match.arg(method)
  counts <- switch(method,
         edgeR = normByDETools(counts, FUN=edgeRnormHelper),
         DESeq2 = normByDETools(counts, FUN=DESeq2normHelper),
         RUVs = normByRUVs(counts, ...),
         fpkm = getFPKM(counts, ...),
         ashr = )
  counts
}

normByDETools <- function(counts, FUN){
  stopifnot(is.function(FUN))
  if("RPFs" %in% names(counts)){
    counts$RPFsRawCounts <- counts$RPFs
    RPFs <- counts$RPFs
    counts$RPFs <- FUN(RPFs)
  }
  
  if("mRNA" %in% names(counts)){
    counts$mRNARawCounts <- counts$mRNA
    mRNA <- counts$mRNA
    counts$mRNA <- FUN(mRNA)
  }
  counts
}

edgeRnormHelper <- function(x){
  if (!requireNamespace("edgeR", quietly=TRUE)) {
    stop("method='edgeR' requires installing the CRAN package 'edgeR'")
  }
  y <- edgeR::DGEList(counts=x)
  y <- edgeR::normLibSizes(y)
  cpm <- edgeR::cpm(y, log = TRUE)
  return(cpm)
}

DESeq2normHelper <- function(x){
  if (!requireNamespace("DESeq2", quietly=TRUE)) {
    stop("method='DESeq2' requires installing the CRAN package 'DESeq2'")
  }
  cd <- data.frame(sample=factor(colnames(x)))
  design <- as.formula('~sample')
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=x,
                                        colData = cd,
                                        design = design)
  dds <- DESeq2::estimateSizeFactors(dds)
  cnt <- DESeq2::counts(dds, normalized=TRUE)
  return(cnt)
}
