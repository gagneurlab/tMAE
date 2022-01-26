merge_scores <- function(data, scores, 
                         populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
                         max_af_cutoff = 0.001,
                         ...){
  

    res <- cbind(data, scores) %>% as.data.table()
    
    # Compute the MAX_AF based on all provided population columns
    # return -1 if only NAs are present (to avoid a warning)
    res$MAX_AF <- apply(res[, ..populations], 1, 
                        FUN=function(x){ max(x, -1, na.rm=TRUE) })
    
    # Replace Inf/-1 with NA
    res[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]
    res[, rare := (MAX_AF <= max_af_cutoff | is.na(MAX_AF))]
    
    return(res)
}

score_data <- function(object, 
    genome_assembly = 'hg19',
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...){
  if("gene_assembly" %in% names(list(...))){
    warning("'gene_assembly' is deprecated. Please use 'genome_assembly' instead.")
    genome_assembly <- list(...)[['gene_assembly']]
  }
  
  if(genome_assembly %in% BiocManager::available("MafDb")){
    mafdb <- .get_mafdb(genome_assembly)
  } else {
    mafdb <- .get_mafdb(switch(genome_assembly, 
      hg19   = "MafDb.gnomAD.r2.1.hs37d5",
      hs37d5 = "MafDb.gnomAD.r2.1.hs37d5",
      hg38   = "MafDb.gnomAD.r2.1.GRCh38",
      GRCh38 = "MafDb.gnomAD.r2.1.GRCh38",
      stop("Please provide a supported genome assembly version. Not: ", genome_assembly)
    ))
  }
    
  if(!all(populations %in% populations(mafdb))){
    stop("Please provide only populations provided by gnomAD!")
  }
  
  # Add score of all, African, American, East Asian and Non-Finnish European
  pt <- score(mafdb, object, pop = populations) %>% as.data.table()
  colnames(pt) <- populations
  
  return(pt)
}

.get_mafdb <- function(pkg_name){
  if(!requireNamespace(pkg_name, quietly=TRUE)){
    warning("The given MafDb is not installed: '", pkg_name, "'. We will do it now!")
    if(!requireNamespace("BiocManager", quietly=TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install(pkg_name, ask=FALSE)
  }
  
  mafdb <- getFromNamespace(pkg_name, pkg_name)
  mafdb
}

MAE_gnomAD <- function(object,...){
  # create GRanges from a MAE count table
  gr <- GRanges(seqnames = object$contig,
                ranges = IRanges(start=object$position, width=1), 
                strand = '*')
  # score the gr
  scores <- score_data(gr,...)
  # merge the scores with the original object
  res <- merge_scores(object,scores,...)
  return (res)
}

GR_gnomAD <- function(object,...){
    
  # create a data table from the GRanges object
  data <- as.data.table(object)
  # score the original GRanges object 
  scores <- score_data(object,...)
  # merge the data.table created and the resulting scores
  res <- merge_scores(data,scores,...)
  return(res)
}

#' @title Add AF to MAE counts
#' @param object a data.frame containing allelic counts.
#' @param list ... any other parameters such as 'max_af_cutoff', or 'populations' within gnomAD
#' @return a data.frame containing allelic counts and minor allele frequencies
#' @export
setMethod("add_gnomAD_AF", signature = "data.table", MAE_gnomAD)

#' @title Add AF to GRanges
#' @description appending the minor allele frequency to GRanges using gnomAD
#' @param object a GRanges object
#' @param list ... any other parameters such as 'max_af_cutoff', or 'populations' within gnomAD
#' @return a data.frame containing GRanges data as well as the minor allele frequencies
#' @export
setMethod("add_gnomAD_AF", signature = "GRanges", GR_gnomAD)
