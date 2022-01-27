#' @title Add AF to GRanges
#' @rdname add_gnomAD_AF
#' @description appending the minor allele frequency to GRanges using gnomAD
#' @param object either a data.table of allelic counts or a GRanges object
#' @param genome_assembly one of "hg19", "hs37d5", "hg38", "GRCh38"
#'                It can also be any full string of a MafDb provided by
#'                \code{\link[GenomicScores]{availableGScores}}.
#' @param max_af_cutoff cutoff for a variant to be considered rare. Default 0.001
#' @param populations The populations to be annotatated.
#' @param ... Used for backwards compatibility (gene_assembly -> genome_assembly)
#' @return a data.frame containing original data as well as the minor allele frequencies
#' @export
setGeneric("add_gnomAD_AF",function(
    object,
    genome_assembly = 'hg19',
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...) standardGeneric("add_gnomAD_AF"))
