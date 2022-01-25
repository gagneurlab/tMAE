test_that("End to end", {
    maeFile <- system.file("extdata", "allelic_counts_HG00187.csv", 
                        package="tMAE", mustWork=TRUE)
    grFile <- system.file("extdata", "GR_HG00187.Rds", 
                        package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(maeFile)
    gr <- readRDS(grFile)
    
    # maeRes
    expect_warning({ maeRes <- DESeq4MAE(maeCounts) }, "NaNs produced")
    expect_warning(add_gnomAD_AF(object=maeRes, gene_assembly='MafDb.ExAC.r1.0.hs37d5', pop="AF"),
            "'gene_assembly' is deprecated")

    # GR
    expect_warning(add_gnomAD_AF(object=gr, gene_assembly='MafDb.ExAC.r1.0.hs37d5', pop="AF"),
            "'gene_assembly' is deprecated")
})
