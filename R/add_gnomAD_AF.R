merge_scores <- function(data, scores, 
                         max_af_cutoff = 0.001,
                         populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
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
    genome_assembly <- switch(genome_assembly, 
      hg19   = "MafDb.gnomAD.r2.1.hs37d5",
      hs37d5 = "MafDb.gnomAD.r2.1.hs37d5",
      hg38   = "MafDb.gnomAD.r2.1.GRCh38",
      GRCh38 = "MafDb.gnomAD.r2.1.GRCh38",
      stop("Please provide a supported genome assembly version. Not: ", genome_assembly) )
    mafdb <- .get_mafdb(genome_assembly)
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

#' @rdname add_gnomAD_AF
#'
#' @examples
#' file <- system.file("extdata", "allelic_counts_HG00187.csv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- fread(file)
#' maeRes <- DESeq4MAE(maeCounts)
#' genome_assembly <- 'hg19'
#' res <- add_gnomAD_AF(maeRes, genome_assembly = genome_assembly, pop="AF")
#'
setMethod("add_gnomAD_AF", signature = "data.table", 
function(
    object, 
    genome_assembly = 'hg19',
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ... ){
  # create GRanges from a MAE count table
  gr <- GRanges(seqnames = object$contig,
                ranges = IRanges(start=object$position, width=1), 
                strand = '*')
  # score the gr
  scores <- score_data(gr,genome_assembly=genome_assembly,populations = populations, max_af_cutoff = max_af_cutoff,...)
  # merge the scores with the original object
  res <- merge_scores(object,scores,populations = populations, max_af_cutoff = max_af_cutoff)
  return (res)
})

#' @rdname add_gnomAD_AF
#'
#' @examples
#' file <- system.file("extdata", "GR_HG00187.Rds", package = "tMAE", mustWork = TRUE)
#' gr <- readRDS(file)
#' genome_assembly <- 'hg19'
#' res <- add_gnomAD_AF(gr, genome_assembly = genome_assembly, pop="AF")
#'
#' 
setMethod("add_gnomAD_AF", signature = "GRanges",
function(
    object, 
    genome_assembly = 'hg19',
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...) {

  # create a data table from the GRanges object
  data <- as.data.table(object)
  # score the original GRanges object 
  scores <- score_data(object,
    genome_assembly = genome_assembly,
    max_af_cutoff = max_af_cutoff,
    populations = populations,
    ...)
  # merge the data.table created and the resulting scores
  res <- merge_scores(data, scores, max_af_cutoff = max_af_cutoff, populations = populations)
  return(res)
})
