#' @title add gnomAD AF to object
#' @description appending the minor allele frequency to GRanges using gnomAD
#' @param object data.table or GRanges object
#' @param ...  Further arguments such as populations, max_af_cutoff, or genome_assembly
#' @return A data.frame containing the added gnomAD MAF
#' @export
setGeneric("add_gnomAD_AF",function(object, ...) standardGeneric("add_gnomAD_AF"))
