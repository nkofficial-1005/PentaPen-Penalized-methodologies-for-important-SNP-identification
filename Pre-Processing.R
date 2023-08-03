##############
# Clear environment
rm(list=ls())
gc()
options(warn=-1)

##############
# Import libraries
library(dplyr)
library(trio)
library(logr)
library(matrixStats)
library(Matrix)
library(dendextend)
library(caTools)
library(caret)
library(tidyverse)
library(grpreg)
##############
# Import data
Geno <- read.pedfile("genotype.ped") 
char.pheno <- read.table("phenotypes.pheno", header = TRUE, stringsAsFactors = FALSE, sep = " ")

##############
# Data Pre-processing
# Re-code the data in ped file
Geno[Geno == 'A'] <- 0  # Converting A to 0
Geno[Geno == 'T'] <- 1  # Converting T to 1
Geno[Geno == 'G'] <- 2  # Converting G to 2
Geno[Geno == 'C'] <- 3  # Converting C to 3

#Convert the phenotype to matrix
y <- matrix(char.pheno$Anthocyanin_22) #Change the phenotype accordingly
rownames(y) <- char.pheno$IID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]

#Imputing SNP null values
for (j in 1:ncol(Geno)) {
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), Geno[, 
                                                                            j])
}
Geno_y <- Geno[index, ] 
pheno_final <- data.frame(famid = rownames(y), y = y)

df <- merge(Geno_y, pheno_final, by = 'famid')
df_final <- df[, 7:214058] #select the data set consisting of SNPs only
df_final <- sapply(df_final, as.numeric)
df_final <- df_final[sample(nrow(df_final)),]
n <- dim(df_final)[[1]] 
d <- dim(df_final)[[2]] 

save(list = ls(),file = "PreprocessedData.RData")