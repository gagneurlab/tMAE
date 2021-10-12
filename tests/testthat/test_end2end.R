test_that("End to end", {
    file <- system.file("extdata", "allelic_counts_HG00187.csv", 
            package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(file)
    expect_warning({ maeRes <- DESeq4MAE(maeCounts) }, "NaNs produced")
    
    if(!requireNamespace("MafDb.gnomAD.r2.1.GRCh38", quietly=TRUE)){
        expect_error(
                add_gnomAD_AF(data=maeRes, genome_assembly='hg38', pop="AF"), 
                "Could not load gnomAD MafDb")
        expect_error(
                add_gnomAD_AF(data=maeRes, genome_assembly='GRCh38', pop="AF"),
                "Could not load gnomAD MafDb")
    } else {
        res <- add_gnomAD_AF(data=maeRes, genome_assembly='hg38', pop="AF")
        res <- add_gnomAD_AF(data=maeRes, genome_assembly='GRCh38', pop="AF")
    }
    
    if(!requireNamespace("MafDb.gnomAD.r2.1.hs37d5", quietly=TRUE)){
        if (!requireNamespace("BiocManager", quietly=TRUE))
            install.packages("BiocManager")
        BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
    }
    
    res <- add_gnomAD_AF(data=maeRes, genome_assembly='hg19', pop="AF")
    res <- add_gnomAD_AF(data=maeRes, genome_assembly='hs37d5', pop="AF")

    expect_is(plotMA4MAE(res), "ggplot")
    expect_is(plotAllelicCounts(res), "ggplot")
})
