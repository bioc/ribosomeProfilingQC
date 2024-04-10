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
#' norm3 <- normBy(cnts, 'vsn')
normBy <- function(counts, method = c('edgeR', 'DESeq2', 'RUVs', 'fpkm', 'vsn'),
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
         vsn = normByDETools(counts, FUN=vsnNormHelper))
  counts
}

normByDETools <- function(counts, FUN){
  stopifnot(is.function(FUN))
  if(all(c("RPFs", 'mRNA') %in% names(counts))){
    if("RPFsRawCounts" %in% names(counts)){
      counts$RPFs <- counts$RPFsRawCounts
    }else{
      counts$RPFsRawCounts <- counts$RPFs
    }
    if("mRNARawCounts" %in% names(counts)){
      counts$mRNA <- counts$mRNARawCounts
    }else{
      counts$mRNARawCounts <- counts$mRNA
    }
    x <- cbind(counts$RPFs, counts$mRNA)
    colnames(x) <- rep(c('RPFs', 'mRNA'), each=4)
    x <- FUN(x)
    RPFs <- x[, seq.int(ncol(counts$RPFs))]
    colnames(RPFs) <- colnames(counts$RPFs)
    counts$RPFs <- RPFs
    mRNA <- x[, -seq.int(ncol(counts$RPFs))]
    colnames(mRNA) <- colnames(counts$mRNA)
    counts$mRNA <- mRNA
    return(counts)
  }
  if("RPFsRawCounts" %in% names(counts)){
    counts$RPFs <- FUN(counts$RPFsRawCounts)
  }else{
    if("mRNA" %in% names(counts)){
      counts$mRNARawCounts <- counts$mRNA
      counts$mRNA <- FUN(counts$mRNA)
    }
  }
  
  if("RPFsRawCounts" %in% names(counts)){
    counts$RPFs <- FUN(counts$RPFsRawCounts)
  }else{
    if("mRNA" %in% names(counts)){
      counts$mRNARawCounts <- counts$mRNA
      counts$mRNA <- FUN(counts$mRNA)
    }
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

vsnNormHelper <- function(x){
  if (!requireNamespace("vsn", quietly=TRUE)) {
    stop("method='vsn' requires installing the CRAN package 'vsn'")
  }
  suppressMessages(Biobase::exprs(vsn::vsn2(x)))
}
