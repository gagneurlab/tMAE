test_that("End to end", {
    file <- system.file("extdata", "allelic_counts_HG00187.csv", 
            package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(file)
    expect_warning({ maeRes <- DESeq4MAE(maeCounts) }, "NaNs produced")
    
    # only run hg38 if existing we do not want to force a download
    pkg_name <- "MafDb.gnomAD.r2.1.GRCh38"
    if(requireNamespace(pkg_name, quietly=TRUE)){
        res <- add_gnomAD_AF(data=maeRes, genome_assembly='hg38', pop="AF")
        res <- add_gnomAD_AF(data=maeRes, genome_assembly='GRCh38', pop="AF")
    }
    
    # expect hg19 to work or to download it on the fly
    res <- add_gnomAD_AF(data=maeRes, genome_assembly='hg19', pop="AF")
    res <- add_gnomAD_AF(data=maeRes, genome_assembly='hs37d5', pop="AF")

    expect_is(plotMA4MAE(res), "ggplot")
    expect_is(plotAllelicCounts(res), "ggplot")
})
