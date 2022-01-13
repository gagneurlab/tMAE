# tMAE
Package containing functions to:

- perform a negative binomial test on allele-specific counts
- add gnomAD minor allele frequencies
- MAplot (FC vs total counts) of allele-specific counts and results
- allelic counts (ALT vs REF)

The results table contains the following information for each SNV:
- refCount: allelic counts of the reference allele
- altCount: allelic counts of the reference allele
- totalCount: total allelic counts
- pvalue: p-value obtained from the negative binomial test using a fix dispersion = 0.05
- padj: Hochberg adjusted p-value
- log2FC: log2 fold change from the negative binomial test (roughly log2(altCount/refCount))
- altRatio: ratio of counts of the alternative allele wrt to the total (altCount / totalCount)
- allele frequencies (AF) of different populations (and the maximum AF), if added
- rare: T/F column if the variant is rare or not, if added
- MAE: T/F if the variant is mono-allelically expressed: (altRatio > 0.8 or altRatio < 0.2) & padj < 0.05
- MAE_ALT: T/F if the alternative allele is expressed: altRatio > 0.8 & padj < 0.05
