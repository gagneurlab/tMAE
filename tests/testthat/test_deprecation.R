test_that("End to end", {
    file <- system.file("extdata", "allelic_counts_HG00187.csv", 
                        package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(file)
    
    expect_warning({ maeRes <- DESeq4MAE(maeCounts) }, "NaNs produced")
    expect_warning(add_gnomAD_AF(data=maeRes, gene_assembly='MafDb.ExAC.r1.0.hs37d5', pop="AF"),
            "'gene_assembly' is deprecated")
})
