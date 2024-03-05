# TransferLearning


# Background
Due to ethnic heterogeneity in genetic structure, polygenic score (PGS) constructed in European populations typically possess poor portability in underrepresented non-European populations. We aim here to explore the predictive performance of the PGS for traits by exploiting existing knowledge of genetic similarity obtained from Europeans. In this work, we refer to the primary population (e.g. EUR) as the auxiliary population and the underrepresented population (e.g. AFR/EAS) as the target population

(https://github.com/BIOstatchen/Transfer-Learning) is implemented in R statistical environment.


# Example
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lmm_PXEM.cpp")


## To describe the SNP effect relationships, we specify a priori functions for SNP effects in Asian populations as follows：β=ωb+δ 
G_EUR=read.table("G_EUR.txt",head=T)
y=read.table("y_EUR.txt",head=T)
X=read.table("X_EUR.txt",head=T)
fit_EUR <- lmm_PXEM(y, X=X, G=G_EUR, PXEM=TRUE, maxIter=1000)
beta_EUR<-as.data.frame(fit_EUR$mub)
beta_EUR=as.matrix(beta_EUR) 
G_EUR=as.matrix(G_EUR)
GRS_EUR<-G_EUR*beta_EUR ## Gb denotes the genetic risk score in the EUR model.



G_EAS=read.table("G_EAS.txt",head=T)
y=read.table("y_EAS.txt",head=T)
X=read.table("X_EAS.txt",head=T)
fit_EAS <- lmm_PXEM(y, X=X, G=G_EAS, PXEM=TRUE, maxIter=1000)
beta_EAS<-as.data.frame(fit_EAS$mub)
beta_EAS=as.matrix(beta_EAS)
G_EAS=as.matrix(G_EAS)
GRS_EAS<-G_EAS*beta_EAS

x<-(G_EAS %*% beta_EUR)
n<-dim(beta_EAS)
fit<-glmnet(x=cbind(x,G_EAS),y=y,alpha=0,penalty.factor=c(0,rep(1,n)),family="gaussian",standardize=TRUE)
coeff<-as.matrix(coef(fit,s="lambda.min")) ##contains ω as a scale parameter and δ as a vector of target-specific influences that follow a normal distribution.
## BETA_EAS=ωb+δ
GRS_EAS_TransferLearning<-G_EAS*BETA_EAS ## denotes the Post transfer Learning genetic risk score in EAS.



//' @param y  response variable
//' @param X  covariates file
//' @param GRS  GRS = G*beta is the GRS information
//' @param G  genotype matrix for SNPs
//' @param maxIter  maximum iteration (default is 1000)

```

# Example1
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("PQLseqMmatglmmAI.cpp")
source("PQLseqMmatglmm.R")


## To describe the SNP effect relationships, we specify a priori functions for SNP effects in Asian populations as follows：β=ωb+δ 
G_EUR=read.table("G_EUR.txt",head=T)
n<-dim(G_EUR)[1]
y=read.table("y_EUR_disease.txt",head=T)
X=read.table("X_EUR.txt",head=T)
k1=(G_EUR %*% t(G_EUR))/n
fit_EUR <- glmmPQL(y=y, X=X, kmat=k1, method="AI",emat=FALSE,model="BMM",verbose=TRUE,maxiter=1000)
u = fit_EUR$lin-(cbind(1,X)%*%fit_EUR$coef)
beta_EUR = solve(t(G_EUR)%*%G_EUR)%*%t(G_EUR)%*%u
beta_EUR=as.matrix(beta_EUR) 
G_EUR=as.matrix(G_EUR)
GRS_EUR<-G_EUR*beta_EUR ## Gb denotes the genetic risk score in the EUR model.



G_EAS=read.table("G_EAS.txt",head=T)
n1<-dim(G_EAS)[1]
y=read.table("y_EAS_disease.txt",head=T)
X=read.table("X_EAS.txt",head=T)
k2=(G_EAS %*% t(G_EAS))/n
fit_EAS <- glmmPQL(y=y, X=X, kmat=k2, method="AI",emat=FALSE,model="BMM",verbose=TRUE,maxiter=1000)
u = fit_EAS$lin-(cbind(1,X)%*%fit_EAS$coef)
beta_EAS = solve(t(G_EAS)%*%G_EAS)%*%t(G_EAS)%*%u
beta_EAS=as.matrix(beta_EAS)
G_EAS=as.matrix(G_EAS)
GRS_EAS<-G_EAS*beta_EAS

x<-(G_EAS %*% beta_EUR)
n<-dim(beta_EAS)
fit<-glmnet(x=cbind(x,G_EAS),y=y,alpha=0,penalty.factor=c(0,rep(1,n)),family="gaussian",standardize=TRUE)
coeff<-as.matrix(coef(fit,s="lambda.min")) ##contains ω as a scale parameter and δ as a vector of target-specific influences that follow a normal distribution.
## BETA_EAS=ωb+δ
GRS_EAS_TransferLearning<-G_EAS*BETA_EAS ## denotes the Post transfer Learning genetic risk score in EAS.



//' @param y  response variable -----(Binary outcomes)
//' @param X  covariates file
//' @param GRS  GRS = G*beta is the GRS information
//' @param G  genotype matrix for SNPs
//' @param maxIter  maximum iteration (default is 1000)

```
# Cite
Wenying Chen, Shuiping Huang and Ping Zeng (2023). Transfer learning prediction of early exposures and genetic risk score on adult obesity in two minority UK Biobank cohorts.$$$#

Wenying Chen, Shuiping Huang and Ping Zeng (2023). Accurate prediction models for the target population through transfer learning shared information from auxiliary populations.$$$#


# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.


# Update
