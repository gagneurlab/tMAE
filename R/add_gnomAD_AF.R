score_data <- function(gr, 
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
  gr_scores <- gscores(mafdb, gr, pop = populations)
  return(gr_scores)
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
  gr_scores <- score_data(gr,genome_assembly=genome_assembly,populations = populations, max_af_cutoff = max_af_cutoff,...)
  
  # merge the population scores with the original object
  scores <- as.data.table(gr_scores)[,..populations]
  res <- cbind(object, scores) %>% as.data.table()
  
  # Compute the MAX_AF based on all provided population columns
  # return -1 if only NAs are present (to avoid a warning)
  res$MAX_AF <- apply(res[, ..populations], 1, 
                      FUN=function(x){ max(x, -1, na.rm=TRUE) })
  
  # Replace Inf/-1 with NA
  res[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]
  res[, rare := (MAX_AF <= max_af_cutoff | is.na(MAX_AF))]
  return (res)
})

#' @rdname add_gnomAD_AF
#'
#' @examples
#' BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
#' file <- system.file("extdata", "GR_HG00187.Rds", package = "tMAE", mustWork = TRUE)
#' gr <- readRDS(file)
#' genome_assembly <- 'hg19'
#' res <- add_gnomAD_AF(gr, genome_assembly = genome_assembly, pop="AF")
#' 
setMethod("add_gnomAD_AF", signature = "GRanges",
function(
    object, 
    genome_assembly = 'hg19',
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...) {

  # score the original GRanges object 
  scores <- score_data(object,
    genome_assembly = genome_assembly,
    max_af_cutoff = max_af_cutoff,
    populations = populations,
    ...)
  return(scores)
})
