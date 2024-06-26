#' Translational Efficiency
#' @description Calculate Translational Efficiency (TE). TE is defined as
#' the ratios of the absolute level of ribosome occupancy divided by RNA levels
#' for transcripts.
#' @param x Output of \link{getFPKM} or \link{normByRUVs}.
#' if window is set, it must be output of \link{coverageDepth}.
#' @param window numeric(1). window size for maximal counts.
#' @param RPFsampleOrder,mRNAsampleOrder Sample order of RPFs and mRNAs.
#' The parameters are used to make sure that the order of RPFs and mRNAs in
#' cvgs is corresponding samples.
#' @param pseudocount The number will be add to sum of reads count to avoid X/0.
#' @param log2 Do log2 transform or not.
#' @param normByLibSize Normalization by library size or not.
#' If window size is provided and normByLibSize is set to TRUE,
#' the coverage will be normalized by library size.
#' @param shrink Shrink the TE or not.
#' @param ... Parameters will be passed to \code{ash} function from \code{ashr}.
#' @return A list with RPFs, mRNA levels and TE as a matrix with
#' translational efficiency
#' @importFrom IRanges RleList IRanges IRangesList viewSums slidingWindows
#' @export
#' @examples
#' \dontrun{
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' RNAs <- dir(path, "mRNA.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' cnts <- countReads(RPFs, RNAs, gtf, level="gene")
#' fpkm <- getFPKM(cnts)
#' te <- translationalEfficiency(fpkm)
#' }
translationalEfficiency <- function(x, window,
                                    RPFsampleOrder, mRNAsampleOrder,
                                    pseudocount=1, log2=FALSE,
                                    normByLibSize=FALSE,
                                    shrink=FALSE,
                                    ...){
  if(!is.list(x)){
    stop("x must be output of countReads, getFPKM, normByRUVs or coverageDepth.")
  }
  if(!all(c("RPFs", "mRNA") %in% names(x))){
    stop("x must be output of countReads, getFPKM, normByRUVs or
         coverageDepth and must contain RPFs and mRNA.")
  }
  RPFs <- x[["RPFs"]]
  mRNA <- x[["mRNA"]]
  shrinkTE <- function(x, ...){
    RPFs <- if(!is.null(x[["RPFsRawCounts"]])){
      x[["RPFsRawCounts"]]
    }else{
      x[["RPFs"]]
    }
    mRNA <- if(!is.null(x[["mRNARawCounts"]])){
      x[["mRNARawCounts"]]
    }else{
      x[["mRNA"]]
    }
    
    if(!all(RPFs==round(RPFs)) &&
       !all(mRNA==round(mRNA))){
      stop('raw counts is missing.')
    }
    if (!requireNamespace("ashr", quietly=TRUE)) {
      stop("type='ashr' requires installing the CRAN package 'ashr'")
    }
    
    if (!requireNamespace("DESeq2", quietly=TRUE)) {
      stop("method='DESeq2' requires installing the CRAN package 'DESeq2'")
    }
    cd <- data.frame(condition=factor(c(rep('RPFs', ncol(RPFs)), 
                                        rep('mRNA', ncol(mRNA)))))
    design <- as.formula('~condition')
    countData <- cbind(RPFs,
                       mRNA)
    mode(countData) <- 'integer'
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=countData,
                                          colData = cd,
                                          design = design)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast=c('condition', 'RPFs', 'mRNA'))
    sebetahat <- res$lfcSE
    betahat <- log2(DESeq2::counts(dds, normalized=TRUE) + pseudocount)
    betahat <- betahat[, seq.int(ncol(RPFs)), drop=FALSE] -
      betahat[, -seq.int(ncol(RPFs)), drop=FALSE]
    
    message("using 'ashr' for TE shrinkage. ",
            "If used in published research, ",
            "please cite: ",
            " Stephens, M. (2016) False discovery rates: ",
            "a new deal. Biostatistics, 18:2. ",
            "https://doi.org/10.1093/biostatistics/kxw041")
    stopifnot('TE is not a matrix with 2 dimension. Please report this bug.'=
                length(dim(betahat))==2)
    fit <- apply(betahat, 2, function(.e){
      ashr::ash(.e, sebetahat,
                mixcompdist="normal", method="shrink", ...)},
      simplify = FALSE)
    TE <- lapply(fit, function(.e) .e$result$PosteriorMean)
    x[['TE']] <- do.call(cbind, TE)
    if(!log2) x[['TE']] <- 2^x[['TE']] - pseudocount
    x
  }
  if(missing(window)){
    if(length(dim(RPFs))!=2 | length(dim(mRNA))!=2){
      stop("x must be output of getFPKM or normByRUVs and
           must contain RPFs and mRNA.")
    }
    if(missing(RPFsampleOrder))
      RPFsampleOrder <- seq.int(ncol(x[["RPFs"]]))
    if(missing(mRNAsampleOrder))
      mRNAsampleOrder <- seq.int(ncol(x[["mRNA"]]))
    id <- intersect(rownames(RPFs), rownames(mRNA))
    if(length(id)==0){return(NULL)}
    x[["RPFs"]] <- RPFs[id, RPFsampleOrder, drop=FALSE]
    x[["mRNA"]] <- mRNA[id, mRNAsampleOrder, drop=FALSE]
    if(shrink){
      x <- shrinkTE(x, ...)
    }else{
      if(normByLibSize){
        normByLib <- function(x){
          xLibSize <- colSums(x, na.rm=TRUE)
          xLibFactor <- mean(xLibSize)/xLibSize
          x <- t(t(x)*xLibFactor)
        }
        RPFs_mRNA <- normByLib(cbind(x[["RPFs"]], x[["mRNA"]]))
        x[["RPFs"]] <- RPFs_mRNA[, seq.int(ncol(x[["RPFs"]])), drop=FALSE]
        x[["mRNA"]] <- RPFs_mRNA[, seq.int(ncol(RPFs_mRNA))[
          -seq.int(ncol(x[["RPFs"]]))], drop=FALSE]
        rm(RPFs_mRNA)
      }
      if(log2){
        x[["TE"]] <- log2(x[["RPFs"]]+pseudocount) - log2(x[["mRNA"]]+pseudocount)
      }else{
        x[["TE"]] <- (x[["RPFs"]]+pseudocount)/(x[["mRNA"]]+pseudocount)
      }
    }
  }else{
    if(!is(RPFs, "cvgd") | !is(mRNA, "cvgd")){
      stop("x must be output of coverageDepth and must contain RPFs and mRNA.")
    }
    if(missing(RPFsampleOrder))
      RPFsampleOrder <- seq.int(length(x[["RPFs"]][["coverage"]]))
    if(missing(mRNAsampleOrder))
      mRNAsampleOrder <- seq.int(length(x[["mRNA"]][["coverage"]]))
    RPFs <- RPFs[["coverage"]][RPFsampleOrder]
    mRNA <- mRNA[["coverage"]][mRNAsampleOrder]
    if(length(mRNA)!=length(RPFs)){
      stop("The length of sample of mRNA is not identical to
           the length of RPFs.")
    }
    if(normByLibSize){
      normByLib <- function(x){
        xLibSize <- lapply(x, sum, na.rm=TRUE)
        xLibSize <- vapply(xLibSize, sum, FUN.VALUE = 0.0, na.rm=TRUE)
        xLibFactor <- mean(xLibSize)/xLibSize
        x <- mapply(x, xLibFactor, FUN=`*`, SIMPLIFY = FALSE)
      }
      RPFs_mRNA <- normByLib(c(RPFs, mRNA))
      RPFs <- RPFs_mRNA[seq_along(RPFs)]
      mRNA <- RPFs_mRNA[seq_along(RPFs_mRNA)[-seq_along(RPFs)]]
      rm(RPFs_mRNA)
    }
    window <- window[1]
    if(round(window)!=window | window < 3){
      stop("window must be a integer no less than 3.")
    }
    ## check all feature_id
    features <- lapply(RPFs, names)
    stopifnot(all(lengths(features)==length(features[[1]])))
    for(i in seq_along(features)){
      stopifnot(all(features[[i]]==features[[1]]))
    }
    features <- features[[1]]

    features2 <- lapply(mRNA, names)
    stopifnot(all(lengths(features2)==length(features2[[1]])))
    for(i in seq_along(features2)){
      stopifnot(all(features2[[i]]==features2[[1]]))
    }
    features2 <- features2[[1]]
    stopifnot(all(features==features2))
    rm(features2)

    cvg <- mapply(RPFs, mRNA, FUN = function(a, b){
      la <- lengths(a)
      lb <- lengths(b)
      stopifnot(identical(la, lb))
      ir <- IRanges(1, la, names = names(la))
      ir <- slidingWindows(ir, window, step = 1L)
      vw.a <- Views(a, ir[names(a)])
      vw.b <- Views(b[names(a)], ir[names(a)])
      sum.a <- viewSums(vw.a, na.rm = TRUE)
      sum.b <- viewSums(vw.b, na.rm = TRUE)
      ratios <- mapply(sum.a, sum.b, FUN = function(r, m){
        if(log2){
          log2(r+pseudocount) - log2(m+pseudocount)
        }else{
          (r+pseudocount)/(m+pseudocount)
        }
      }, SIMPLIFY = FALSE)
      ids <- unlist(lapply(ratios, which.max))
      ratios <- mapply(ratios, ids, FUN=function(value, key){
        value[key]
      })
      rpf <- mapply(sum.a, ids, FUN=function(value, key){
        value[key]
      })
      mrna <- mapply(sum.b, ids, FUN=function(value, key){
        value[key]
      })
      list(ratios=ratios, RPFs=rpf, mRNA=mrna)
    }, SIMPLIFY = FALSE)
    x[["RPFs"]] <- do.call(cbind, lapply(cvg, `[[`, i="RPFs"))
    x[["mRNA"]] <- do.call(cbind, lapply(cvg, `[[`, i="mRNA"))
    x[["TE"]] <- do.call(cbind, lapply(cvg, `[[`, i="ratios"))
  }
  return(x)
}
