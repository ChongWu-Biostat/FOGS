# FOGS (The new version is at [ChongWuLab](https://github.com/ChongWuLab/FOGS))

**Please check our updated version, whic is much easier to run!**

Transcriptome-wide association studies (TWAS) have been recently applied to successfully identify many novel genes associated with complex traits. While appealing, TWAS tend to identify multiple significant genes per locus and many of them may not be causal due to confounding through linkage disequilibrium (LD) among SNPs. Here we introduce a powerful fine-mapping method called FOGS that prioritizes putative causal genes by accounting for local LD. Check the following manuscript for details:

Chong Wu and Wei Pan. A powerful fine-mapping method for transcriptome-wide association studies. *Under review*.



## Installation

Before applying FOGS, we need to prepare the supporting data and install some necessary packages.

1. Download and unpack the LD reference panel. For illustration, we can use the following command to download 1000 Genomes reference:

   > ```
   > wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
   > tar xjvf LDREF.tar.bz2
   > ```

2. Download and prepare the weights from the TWAS(http://gusevlab.org/projects/fusion/) or PrediXcan (http://predictdb.org/) website

3. Download loci definition from LDetect (https://bitbucket.org/nygcresearch/ldetect/src)

4. Launch R and install required libraries:

> ```R
> install.packages('optparse','data.table','matlib','Rcpp','RcppArmadillo','bigmemory','mvtnorm','MASS','dplyr','GenomicRanges','mvnfast')
> install.packages('plink2R-master/plink2R/',repos=NULL)
> if (!require("devtools"))
> install.packages("devtools")
> devtools::install_github("ChongWu-Biostat/aSPU2")
> devtools::install_github("gabraham/plink2R/plink2R")
> ```

5. Install plink (https://www.cog-genomics.org/plink/1.9/) and add the path to the global path by: export PATH=$PATH:/plink_directory. Note that you can revise the .bashrc to add the path prenatally, but we recommend add it every time to avoid some software conflicting.



## Main function

Our FOGS method involves two steps. To help researcher develop some new methods, we summarize our core functions as follows:

- GWAS summary data based ridge regression:  JointRidge.R

  > ```R
  > source("jointRidge.R")
  > 
  > B <- sumstat.tmp2[, "beta"] #effect size
  > S <- sumstat.tmp2[, "se"] #corresponding standard deviation
  > N <- sumstat.tmp2[, "N"] #corresponding smaple size
  > 
  > cur.genos <- cur.genos - matrix(rep(colMeans(cur.genos), each = dim(cur.genos)[1]), dim(cur.genos)[1], dim(cur.genos)[2])
  > XX <- cov(cur.genos) # reference panel
  > 
  > lambda <- 0.1
  > res <- JointRidge(B, S, N, XX, lambda)
  > beta <- res$beta #joint effect size
  > se <- sqrt(diag(res$cov)) # corresponding standard deviation
  > cond.Z <- beta / se # conditional z score
  > ```

- We calculate the conditional z score for each eQTL derived SNP, and then combine conditional Z score by aSPU test:

  > calc.pvalue(cond.Z) # cur.Z conditional Z score

## Pipeline

To help researchers better use the software and conduct the fine-mapping with other TWAS-related methods. We share our pipeline for our real data analysis.

- [ ] Step 1: prepare the eQTL-derived weights. As FOCUS,  we used the causal tissue-based eQTL weights when available; otherwise, we included the eQTL weights with the best accuracy across all other tissues (pre_process_weight_elasticnet.R).
- [ ] Step 2: prepare the pruned-SNP set information (pre_process_prepare_SNP_inf.R)
- [ ] Step 3: conduct FOGS (FOGS_scz.R)

