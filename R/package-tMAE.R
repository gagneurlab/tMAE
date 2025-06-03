#'
#' This is the import mapping for the tMAE package
#' 
#' @noRd
#' 
#' @name tMAE
#'
#' @import data.table
#' 
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot geom_point theme_bw scale_y_log10 scale_x_log10 labs scale_color_manual theme aes element_blank geom_abline
#' @importFrom utils getFromNamespace install.packages
#' 
#' @importFrom BiocGenerics estimateSizeFactors plotDispEsts
#' 
#' @importFrom DESeq2 normalizationFactors normalizationFactors<-
#'          sizeFactors sizeFactors<- counts counts<-
#'          DESeqDataSetFromMatrix DESeqDataSet
#'          makeExampleDESeqDataSet show fpkm fpm
#'          estimateSizeFactorsForMatrix replaceOutliers
#'          dispersions dispersions<- nbinomWaldTest results
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicScores populations seqnames gscores
#' @importFrom BiocGenerics score start pos
#' @importFrom IRanges IRanges
#' 
NULL

globalVariables(c(
        ".",
        "..populations",
        "altCount",
        "altRatio",
        "FC",
        "hgncid",
        "MAX_AF",
        "padj",
        "rare",
        "refCount",
    	"Significant",
        "totalCount"),
    package="tMAE")


