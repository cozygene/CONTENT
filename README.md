
# CONTENT 
### A method for multi-context genetic modeling of transcription regulation

Currently included in the repository is the software for running CONTENT as well as the context-by-context (classic TWAS/predXcan) approach. The implementation of both methods takes advantage of fast implementations of penalized linear models through package `bigstatsr`[1]. Future plans for development include a pipeline and instructions for conducting your own TWAS with hierarchical FDR[2,3]. In the meantime, feel free to reach out to us at mjthompson (at) ucla (dot) edu if you would like us to help you with your own TWAS using CONTENT and/or the context-by-context approach as well as on GTEx or our single-cell PBMC RNAseq dataset.

## Running CONTENT and the context-by-context approach

CONTENT is run as an Rscript via command line, one run per gene (and all its measured contexts). This script will generate weights and cross-validated predictors for both the context-by-context approach and CONTENT, assuming there are at least 3 contexts for the CONTENT approach to be identifiable. The required arguments are as follows:

 - `--X_file` A genotype matrix with the set of individuals in all Y_files and their corresponding cis-SNPs. Row names are individuals, no column names. 
 - `--Y_file_dir` Directory containing only the expression files for each context. The individuals (row names) must be a subset of those in the genotype (X_file). Contains an unnamed column with the expression---or if cov_file_dir is NULL, the residual expression for a given context.
 - `--cov_file_dir` Specify this directory if you want to residualize the expression. The covariates are context-specific and should share the same prefix as the expression files such that file $i$ of the expression corresponds to file $i$ of the covariates. E.g. "Liver_expression" and "Liver_covariates". If NULL, CONTENT will run assuming the expression has already been residualized. Needs to (1) be the same row-length as the corresponding expression file (2) contain the same rownames (3) contain no column names.
 - `--out_dir` Where do save the cross-validated predictors and their corresponding performance statistics. Use the same directory for all genes, and it will contain the predictors as well as their pvalues and $R^2$ values for each gene for convenient reporting.
 - `--snps` A tab-delimited txt file containing information from your .bed or other genotype file. Contains 6 columns and the number of rows corresponds to the number of snps. No column names or row names. **Col2 must be rsIDs (or matching your GWAS summary statistic SNP names) in order for FUSION/TWAS to work**
	 1. chromosome
	 2. rsID
	 3. location CM (not exactly necessary for TWAS)
	 4. location/locus on chromosome
	 5. allele1
	 6. allele2
-	`--twas_dir` Where to save the TWAS weights. We recommend using the same directory for all genes and it will conveniently store the results (such that all significantly cross-validated predicted genes will have their weights in the same directory and we can point FUSION/TWAS to this single directory).
-	`--gene_name` Add this prefix to all the saved results files, this is necessary to distinguish your results.

There are also additional arguments that one can change:

 - `--verbose`Show progress and output as the script runs.
 - `--signif` Manually override the default (0.1) nominal prediction pvalue of when a gene should have weights built for use in TWAS. Our analysis shows there remains power to be gained by setting this somewhat liberally at 0.1
 - `--content_alpha` Linear model regularization constant for CONTENT. Default is elastic net (0.5), LASSO is 1 and ridge-like is 0.001.
 - `--tissue_alpha` Linear model regularization constant for the context-by-context approach.
 - `--num_folds` Number of folds for cross-validation
 - --`scale_exp` Center and scale the gene expression

	
[1] Florian Prive, Hugues Aschard, et al. “Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr.” In:Bioinformatics34.16 (2018), pp. 2781–2787.

[2] C B Peterson, M Bogomolov, et al. “TreeQTL: hierarchical error control for eQTL findings.” In:Bioinformatics 32.16 (2016), pp. 2556–2558.

[3] Michael Wainberg, Nasa Sinnott-Armstrong, et al. “Opportunities and challenges for transcriptome-wide association studies”. In:Nature Genetics 51.4 (2019), pp. 592–599.


