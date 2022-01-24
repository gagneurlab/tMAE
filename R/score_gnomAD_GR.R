#' Add allele frequencies from gnomAD
#'
#' @description Helper Method for add_gnomAD_AF.R generates the gnomAD scoring table given a GenomicRanges object.
#' @author Vicente Yepez and Nicholas Smith
#' @param gr GenomicRanges object
#' @param genome_assembly either 'hg19/hs37d5' or 'hg38/GRCh38' indicating the genome assembly of the variants.
#'                It can also be any full string of a MafDb provided by 
#'                \code{\link[GenomicScores]{availableGScores}}.
#' @param max_af_cutoff cutoff for a variant to be considered rare. Default is .001.
#' @param pops The population to be annotated.
#' @return A data.table with the scores of each position in the gr object using the mafdb database
#' @export
#' 

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
