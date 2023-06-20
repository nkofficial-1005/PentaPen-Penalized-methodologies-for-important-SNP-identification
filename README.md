#  PentaPen: Combining Penalized Models for Identification of Important SNPs on Whole-genome Arabidopsis Thaliana Data

#  Table of content

* [Abstract](#abstract)

* [Technologies](#technologies)

* [Data](#data)

* [Step by step implementation](#step-by-step-implementation)

#  Abstract

Genome-Wide Association Study (GWAS) is the discovery of  an association between certain variations in the genetic code (genome) and a certain physical trait (phenotype). Single Nucleotide Polymorphisms (SNPs) are the most abundant form of common simple DNA variants. In bioinformatics studies, one of the most challenging processes to carry out association tests is finding significant SNPs in high-dimensional data. This problem can be potentially solved by feature selection using statistical and machine learning algorithms. This research focuses on proposing a new pipeline that works best to find the most significant SNPs and comparing the penalized techniques in two phases (using all and filtered SNPs) on  Arabidopsis thaliana  data. The SNPs discovered using the proposed pipeline are further compared with the results of GWAS software, GAPIT. The penalized methodologies used  in this study include Ridge, elastic net, LASSO, Group LASSO, and  Sparse Group LASSO (SGL). Results of binary, continuous, and categorical phenotype indicate that when all the SNPs are used as the input to the penalized methodologies, the LASSO and elastic net outperforms Ridge from two experimental data sets. A pipeline is proposed that uses a filtering strategy called SNP Pooling to narrow the number of SNPs. It takes the union of the SNPs output from the methods Ridge, LASSO, and elastic net, and then sends them as the input for Group LASSO and  Sparse Group LASSO (SGL). Results of the proposed pipeline indicate  that it outperforms the single penalized methodologies and the output SNPs are more beneficial for follow-up analysis. The study provides a rigorous comparison of penalized methodologies for interested researchers  about the best available methods for feature selection. It also provides  a seamless pipeline for significant SNP identification while reducing the computational work and selections of the preferred method.

Key Words: Genomic Wide Association Study  ·  Single Nucleotide Polymorphism  ·  Feature Selection  ·  Machine Learning  ·  High Dimensional Data.

#  Technologies

Software: R Version 4.2.2

Operating System: Linux 5.4.0-135-generic x86_64

Cloud Server: TRU Data Science

#  Data


Two  Arabidopsis thaliana  data, AtPolyDB and F1, are used for this study. They are obtained from easygwas websites: https://easygwas.ethz.ch/data/public/dataset/view/1/ and https://easygwas.ethz.ch/data/public/dataset/view/42/. The AtPolyDB dataset has 1307 samples with 214051 SNPs (or features) and the F1 data set has 372 samples with 204753 SNPs. Both data sets contain three files: (a) PED file, (b) PHENO file, and (c) MAP file. The chosen phenotypes had three different data types: (a) Binary, (b) Continuous, and (c) Categorical.

#  Step-by-step implementation

The new pipeline includes below steps,

<b>Input:</b>  Genotype .ped file and Phenotype .pheno file  
1. Re-code the chromosomal nucleotide to numeric values to form a binary marker data followed by creating a design matrix of dimensions n  ×  p.
2. Remove null values from the phenotype data and match them with marker data.  
3. Impute SNPs with null values with the mean across all samples in the marker data set.  
4. For a 5-fold CV, repeat the steps  
a. Split the data into training (80%) and testing (20%) folds.  
b. Use glmnet to train and validate ridge, lasso, and elastic net.
	1. Predict the phenotype value for both training and testing  folds using  λ  within 1 standard error of the minimum obtained by an inner 5-fold CV.  
	2. Record the appropriate performance metric using the optimal cutoff to optimize the metric.  
	3. Record the potentially significant SNPs from each method. The significant SNPs are those with coefficients higher than the mean of the absolute value of the coefficients.

	c. Filter SNPs by taking the union of SNPs from the ridge, LASSO, and elastic net.  
d. Create groups of SNPs using Hierarchical Clustering.  
e. Utilize filtered SNPs to train and validate Group Lasso and SGL using R functions grplasso and SGL.
	1. Predict the phenotype value for both training and testing folds using  λ  within 1 standard error of the minimum obtained by an inner 5-fold CV.  
	2. Record the appropriate performance metric using the optimal cutoff to optimize the metric.  
	3. Take the union of the potentially significant SNPs from both Group Lasso and SGL. The significant SNPs are those with coefficients higher than a cutoff (mean of the absolute value of the coefficients).

<b>Output:</b>  The significant SNPs (union of selected SNPs from Group LASSO and SGL) for each phenotype
