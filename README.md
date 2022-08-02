# Ensemble-Learning-Based-Multimarker-Association-Test
Testing association of a group of markers with continuous and case-control phenotypes

The following R program can be used for testing association of a group of markers with continuous and case-control traits. (For methodological details please refer to: Padhukasahasram B, Reddy CK, Levin AM, Burchard EG, Williams LK (2015) Powerful Tests for Multi-Marker Association Analysis Using Ensemble Learning. PLoS ONE 10(11): e0143489. doi:10.1371/journal.pone.0143489)

https://github.com/pkbadri/ensemble-learning-based-multimarker-association-test/blob/main/multimarker-association-version2.0.R

The following are the command line parameters to be specified in this same order:

-snps Number of markers.

-features Number of markers to use in predictive model. Must be <= snps.

-samples Number of samples.

-trait 1 for continuous traits and 2 for case-control data.

-covs Number of covariates to adjust for.

-snpfile Filename storing SNP genotype data (e.g. 0, 1, or 2). Numbers must be separated by spaces and rows represent samples while columns represent markers.

-phenofile Filename storing phenotype information. Each line represents the phenotype value for samples in the same order as in snpfile.

-covfile Filename storing information about covariates. Each row should have values of covariates for a sample separated by spaces. Samples should be in the same order as in snpfile.

For example:

Rscript multimarker-association-version2.0-faster.R -snps 12 -features 3 -samples 2790 -trait 1 -covs 3 -snpfile snps -phenofile pheno -covfile covariates



















# Ensemble-Learning-Based-Multimarker-Interactions-Test


The following R program can be used to conduct a multi-marker test for interactions.

https://github.com/pkbadri/ensemble-learning-based-multimarker-association-test/blob/main/multimarker-interactions-test.R

Usage is same as before. For example:
Rscript multimarker-interactions-test.R -snps 12 -features 3 -samples 2790 -trait 1 -covs 3 -snpfile snps -phenofile pheno -covfile covariates
