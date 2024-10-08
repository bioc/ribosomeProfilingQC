#' Codon usage bias
#' @description Calculate the codon usage for the reads in the identified CDSs.
#' And then compared to the reference codon usage.
#' @param RPFs Bam file names of RPFs.
#' @param gtf GTF file name for annotation or a TxDb object.
#' @param genome A BSgenome object.
#' @param bestpsite P site postion.
#' @param readsLen Reads length to keep.
#' @param anchor 5end or 3end. Default is 5end.
#' @param ignore.seqlevelsStyle Ignore the sequence name style detection or not.
#' @param summary Return the summary of codon usage bias or full list.
#' @param removeDuplicates Remove the PCR duplicates or not. Default TRUE.
#' @param ... Parameters pass to
#' \link[txdbmaker:makeTxDbFromGFF]{makeTxDbFromGFF}
#' @return A list of data frame of codon count table if summary is TRUE.
#'  list 'reads' means the counts by raw reads.
#'  list 'reference' means the counts by sequence extracted from reference by
#'  the coordinates of mapped reads.
#'  Otherwise, return the counts (reads/reference) table for each reads.
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings trinucleotideFrequency GENETIC_CODE AMINO_ACID_CODE
#' @export
#' @examples
#' path <- system.file("extdata", package="ribosomeProfilingQC")
#' RPFs <- dir(path, "RPF.*?\\.[12].bam$", full.names=TRUE)
#' gtf <- file.path(path, "Danio_rerio.GRCz10.91.chr1.gtf.gz")
#' library(BSgenome.Drerio.UCSC.danRer10)
#' cb <- codonBias(RPFs[c(1,2)], gtf=gtf, genome=Drerio)
codonBias <- function(RPFs, gtf, genome,
                      bestpsite=13,
                      readsLen=c(28,29),
                      anchor="5end",
                      ignore.seqlevelsStyle=FALSE,
                      summary=TRUE,
                      removeDuplicates = TRUE,
                      ...){
  yieldSize <- 10000000
  stopifnot(is.logical(summary))
  summary <- summary[1]
  stopifnot(is.character(gtf)||is(gtf, "TxDb"))
  stopifnot(inherits(genome, c("DNAStringSet", "BSgenome")))
  anchor <- match.arg(anchor, choices = c("5end", "3end"))
  stopifnot(is.character(RPFs))
  stopifnot(is.numeric(readsLen))
  stopifnot(is.numeric(bestpsite))
  txdb <- prepareTxDb(gtf, 'gtf', ...)
  cds <- cdsBy(txdb, by = 'tx', use.names = TRUE)
  refSeq <- getSeq(genome, cds)
  refSeq <- vapply(refSeq, paste, FUN.VALUE=character(1L), collapse='')
  res <- lapply(RPFs, function(f){
    bamfile <- BamFile(file = f, yieldSize = yieldSize)
    pc <- getPsiteCoordinates(
      bamfile, bestpsite=bestpsite,
      anchor = anchor,
      param = ScanBamParam(what=c("qwidth", "seq"),
                           tag=character(0),
                           flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                            isUnmappedQuery = FALSE,
                                            isNotPassingQualityControls = FALSE,
                                            isSupplementaryAlignment = FALSE)))
    if(removeDuplicates){
      pc <- pc[!duplicated(paste(seqnames(pc), start(pc), pc$qwidth))]
    }
    pc <- pc[pc$qwidth %in% readsLen]
    pc <- shiftReadsByFrame(pc, txdb,
                                ignore.seqlevelsStyle=ignore.seqlevelsStyle)
    pc <- pc[!is.na(pc$tx_name)]
    if(anchor=="5end"){
      ## -3 is docking position, smaller than that should not be ready
      pc <- pc[pc$posToStop>=0 & pc$position + pc$Psite>=-3]
    }else{
      ## P site is stop coden, the transcripts is leaving
      pc <- pc[pc$position>=0 & pc$posToStop + 3 >= pc$Psite]
    }
    
    #the sequence in bam file should be the raw reads, no need to do reverseComplement.
    pc$char <- substring(pc$seq, ifelse(pc$Psite>pc$position, pc$Psite-pc$position, 1))
    pcs <- ifelse(pc$Psite>pc$position, 0, pc$position-pc$Psite+1)
    pc$ref <- substr(refSeq[pc$tx_name],
                     pcs+1,  pcs + nchar(pc$char))
    pc$char <- substr(pc$char, 1, nchar(pc$ref))
    
    ## position, relative position from start codon
    ## posToStop, relative position from stop codon
    seqStart <- pc$position %% 3
    seqWidth <- pc$qwidth - seqStart
    seqWidth <- seqWidth - (seqWidth %% 3)
    
    pc$char <- substr(pc$char, seqStart+1, seqStart + seqWidth)
    pc$ref <- substr(pc$ref, seqStart+1, seqStart + seqWidth)

    seqCodonUsage <- trinucleotideFrequency(DNAStringSet(pc$char), step = 3)
    refCodonUsage <- trinucleotideFrequency(DNAStringSet(pc$ref), step = 3)
    if(summary){
      ## remove the PCR duplicates or not?
      seqCodonUsage <- colSums(seqCodonUsage)
      refCodonUsage <- colSums(refCodonUsage)
      return(data.frame(reads=seqCodonUsage, reference=refCodonUsage))
    }else{
      codonUsage <- paste(seqCodonUsage, refCodonUsage, sep=',')
      dim(codonUsage) <- dim(seqCodonUsage)
      mcols(pc) <- cbind(mcols(pc), codonUsage)
      return(pc)
    }
  })
  names(res) <- basename(RPFs)
  if(summary){
    reshapeData <- function(res, column){
      dat <- do.call(cbind, lapply(res, function(.ele) .ele[, column]))
      rn <- rownames(res[[1]]) 
      dat <- cbind(codon = rn,
                   AAcodon = GENETIC_CODE[rn],
                   Abbreviation = AMINO_ACID_CODE[GENETIC_CODE[rn]],
                   data.frame(dat))
      dat
    }
    reads <- reshapeData(res, 'reads')
    reference <- reshapeData(res, 'reference')
    res <- list(reads=reads, reference=reference)
  }
  res
}
