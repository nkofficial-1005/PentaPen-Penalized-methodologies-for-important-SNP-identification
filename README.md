#  PentaPen: Combining Penalized Models for Identification of Important SNPs on Whole-genome Arabidopsis Thaliana Data

#  Table of content

* [Abstract](#abstract)

* [Technologies](#technologies)

* [Data](#data)

* [Step by step implementation](#step-by-step-implementation)
  
* [Instructions to Run Code](#instructions-to-run-code)
  
#  Abstract

Genome-Wide Association Study (GWAS) is the discovery of  an association between certain variations in the genetic code (genome) and a certain physical trait (phenotype). Single Nucleotide Polymorphisms (SNPs) are the most abundant form of common simple DNA variants. In bioinformatics studies, one of the most challenging processes to carry out association tests is finding significant SNPs in high-dimensional data. This problem can be potentially solved by feature selection using statistical and machine learning algorithms. An improved penalized-method-based workflow, PentaPen, is  developed to find important SNPs combined with different penalized models. "Penta" and "Pen" are abbreviated for "five" and "penalized models" respectively. PentaPen is a classifier and regressor which is developed for SNP identification. It aims to minimize the value of the loss function while simultaneously optimizing performance metrics. Firstly, all the SNPs of whole-genome SNP data are utilized for training Ridge, LASSO, and Elastic Net. The union of the output SNPs from these three models is taken as the selected SNPs known as SNP Pooling. Secondly, the selected SNPs are sent to train Group LASSO and SGL, and the union of the output SNPs from the two is the final output of the PentaPen. Finally, an aggregated model is developed by combining the predictions of all five penalized models; this model is used to calculate the performance metrics of PentaPen. The proposed workflow aims to enhance the confidence of the selected SNPs by leveraging the beneficial properties of five penalized methodologies. As a result, combining multiple penalized models can improve performance by reducing over-fitting as compared to using one model. The workflow for SNP identification can also increase the confidence in choosing an SNP set as it is more probable to select informative SNPs than using a single penalized model. Hence, the union of SNPs from Group LASSO and SGL allows for further analysis and selection of a reduced number of SNPs since SGL's sparsity group-wise and within-group results in too few or no SNPs.

Key Words: Genomic Wide Association Study  ·  Single Nucleotide Polymorphism  ·  Feature Selection  ·  Machine Learning  ·  High Dimensional Data.

#  Technologies

Software: R Version 4.2.2 and R Version 3.6.3

Operating Systems: Linux 5.4.0-135-generic x86_64 and Linux 5.4.0-150-generic x86_64

Cloud Servers: TRU Data Science and Compute Canada

#  Data

Two  Arabidopsis thaliana  data, AtPolyDB and F1, are used for this study. They are obtained from easygwas websites: https://easygwas.ethz.ch/data/public/dataset/view/1/ and https://easygwas.ethz.ch/data/public/dataset/view/42/. The AtPolyDB dataset has 1307 samples with 214051 SNPs (or features) and the F1 data set has 372 samples with 204753 SNPs. Both data sets contain three files: (a) PED file, (b) PHENO file, and (c) MAP file. The chosen phenotypes had three different data types: (a) Binary (Anthocyanin), (b) Continuous (Width and DTF), and (c) Categorical (Germination Days).

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

# Instructions to Run Code
1. Download all R files and keep them in the same directory.
2. In _Pre-Processing.R_ file  
a. Input the Genotype .ped file and Phenotype .pheno file.  
b. Now, Select the appropriate phenotype data type in Line 33.
c. Save the list of objects as PreprocessedData.RData. (Code included)
	4. In _PentaPenCode.R_ file
   a. Load all libraries and the PreprocessedData. (Code included)
   b. Now, run the 5-fold CV and parallel computing within the for loop to get the aggregated model (PentaPen).
(The functions to train five penalized models are defined in _BeforeSNPPool.R_ and _AfterSNPPool.R_ files. They are called during parallel computation.)
   d. Finally, run the aggregated model to get the final results.
(Comment out the evaluation metrics according to the chosen phenotype.)
