#' Add allele frequencies from gnomAD
#'
#' @description Add allele frequency information from gnomAD.
#' @author Vicente Yepez
#' @param data A data.frame containing allelic counts.
#' @param genome_assembly either 'hg19/hs37d5' or 'hg38/GRCh38' indicating the genome assembly of the variants.
#'                It can also be any full string of a MafDb provided by 
#'                \code{\link[GenomicScores]{availableGScores}}.
#' @param max_af_cutoff cutoff for a variant to be considered rare. Default is .001.
#' @param pops The population to be annotated.
#' @param ... Used for backwards compatibility (gene_assembly -> genome_assembly)
#' @return A data.table with the original contents plus columns containing allele frequencies from different gnomAD populations.
#' @export
#' 
#' @examples
#' file <- system.file("extdata", "allelic_counts_HG00187.csv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- fread(file)
#' maeRes <- DESeq4MAE(maeCounts)
#' 
#' # define the assembly/MafDb you want e.g. hg19, MafDb.gnomAD.r2.1.hs37d5, or MafDb.ExAC.r1.0.hs37d5
#' \dontrun{
#' genome_assembly <- 'hg19' 
#' maeRes <- add_gnomAD_AF(maeCounts, genome_assembly = genome_assembly, pop="AF")
#' }
#' 

# use score_gnomAD_GR as helper to add gnomAD frequencies to data
add_gnomAD_AF <- function(data,
                       genome_assembly = c('hg19', 'hs37d5', 'hg38', 'GRCh38'),
                       max_af_cutoff = .001,
                       pops = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
                       ...){
  
  if("gene_assembly" %in% names(list(...))){
    warning("'gene_assembly' is deprecated. Please use 'genome_assembly' instead.")
    genome_assembly <- list(...)[['gene_assembly']]
  }
  # Transform data into GRanges object
  gr <- GRanges(seqnames = data$contig,
                ranges = IRanges(start=data$position, width=1), 
                strand = '*')
  
  # add scores to data table
  score_table <- score_gnomAD_GR(gr,genome_assembly,max_af_cutoff,pops)
  res <- cbind(data, score_table) %>% as.data.table()
  
  # Compute the MAX_AF based on all provided population columns
  # return -1 if only NAs are present (to avoid a warning)
  res$MAX_AF <- apply(res[, ..pops], 1, 
                      FUN=function(x){ max(x, -1, na.rm=TRUE) })
  
  # Replace Inf/-1 with NA
  res[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]
  res[, rare := (MAX_AF <= max_af_cutoff | is.na(MAX_AF))]
  
  return(res)
}

# add gnomAD frequencies to granges
score_gnomAD_GR <- function(gr, 
    genome_assembly = c('hg19', 'hs37d5', 'hg38', 'GRCh38'),
    max_af_cutoff = .001,
    pops = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...){
  
  if(genome_assembly %in% BiocManager::available("MafDb")){
    mafdb <- .get_mafdb(genome_assembly)
  } else {
    genome_assembly <- match.arg(genome_assembly)
    mafdb <- .get_mafdb(switch(genome_assembly, 
      hg19   = "MafDb.gnomAD.r2.1.hs37d5",
      hs37d5 = "MafDb.gnomAD.r2.1.hs37d5",
      hg38   = "MafDb.gnomAD.r2.1.GRCh38",
      GRCh38 = "MafDb.gnomAD.r2.1.GRCh38",
      stop("Please provide a supported genome assembly version. Not: ", genome_assembly)
    ))
  }
    
  if(!all(pops %in% populations(mafdb))){
    stop("Please provide only populations provided by gnomAD!")
  }
  
  # Add score of all, African, American, East Asian and Non-Finnish European
  pt <- score(mafdb, gr, pop = pops) %>% as.data.table()
  colnames(pt) <- pops
  
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
