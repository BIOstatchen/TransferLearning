# TransferLearning


# Background
Due to ethnic heterogeneity in genetic structure, genetic risk scores (GRS) constructed in European populations typically possess poor portability in underrepresented non-European populations. We aim here to explore the predictive performance of the GRS for traits by exploiting existing knowledge of genetic similarity obtained from Europeans. In this work, we refer to the primary population (e.g. EUR) as the auxiliary population and the underrepresented population (e.g. AFR/EAS) as the target population

(https://github.com/BIOstatchen/Transfer-Learning) is implemented in R statistical environment.


# Example

library(glmnet)
source("lmm_PXEM_Rcpp.R")
source("estimate_beta.R")
BETA=read.table("BETA.txt",head=T) ## BETA is the effect size of SNPs in base population from summary statistics and SE is their standard error.
SE=read.table("SE.txt",head=T)
R=read.table("R.txt",head=T) ## R is the linkage disequilibrium (LD) matrix which can be calculated with genotypes of population-matched individuals from external reference panels such as the 1000 Genomes Project in the base population. We calculate R in a shrinkage fashion R = 0.95 * as.matrix(cor(G_base_population)) + diag(1-0.95, nrow=nSNPs, ncol=nSNPs)
G=read.table("G.txt",head=T)
y=read.table("y.txt",head=T)
X=read.table("X.txt",head=T)
BETA=as.matrix(BETA)
SE=as.matrix(SE)
R=as.matrix(R)
G=as.matrix(G)
y=as.matrix(y)
X=as.matrix(X)
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)

//' @param y  response variable
//' @param X  covariates
//' @param g  GRS = G*weight is the GRS information
//' @param G  genotype matrix for GWAS
//' @param maxIter  maximum iteration (default is 1000)















# Cite
Wenying Chen, Shuiping Huang and Ping Zeng (2022). Transfer learning prediction of early exposures and genetic risk score on adult obesity in two minority UK Biobank cohorts.$$$#



# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.


# Update
