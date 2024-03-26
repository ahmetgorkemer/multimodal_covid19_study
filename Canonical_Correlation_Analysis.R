# Sample script for canonical correlation analysis of Radiomic Features and Viral Genome Word2Vec Encoding. This script can be extended to analyze the correlation between other pairs.

library(PMA)
library(ggplot2)
library(dplyr)
library("data.table")
library(tidyr)

file <- '...csv'

data <- read.table(file, sep=",",stringsAsFactors = FALSE)
data <- data %>% row_to_names(row_number = 1)     # All data modalities 

data$Age = as.numeric(as.character(data$Age))
data$Charlson_Comorbidity_Score = as.numeric(as.character(data$Charlson_Comorbidity_Score))

c <- as.numeric(23)     # Initial column numbers of relevant data modalities. C: Clinical Data; B: Viral Genome Binary Encoding; R: Radiomics Features; V: Viral Genome Word2Vec Encoding; L: Laboratory Data
b <- as.numeric(57)
r <- as.numeric(629)
v <- as.numeric(2205)
h <- as.numeric(2505)
l <- as.numeric(2519)

x_rv <- data %>% drop_na(r,v,icu)   # An example analysis: Canonical Correlation Analysis of Viral Genome Word2Vec Encoding and Radiomics Features
r_rv <- x_rv[,c(r:(v-1))]
v_rv <- x_rv[,c(v:(h-1))]

x <- apply(as.matrix(r_rv),2,as.numeric)
z <- apply(as.matrix(v_rv),2,as.numeric)

# Filter by variance
X_ind = which(apply(x,2,var) > 0.001) 
Z_ind = which(apply(z,2,var) > 0.001) 
x_fil = x[,X_ind]
z_fil = z[,Z_ind]

set.seed(3)
rv_cca_fit_perm = CCA.permute(x_fil,z_fil,
                              typex="standard",
                              typez="standard",
                              nperms=25,trace=FALSE)

rv_cca_fit = CCA(x_fil,z_fil,typex = 'standard',typez='standard',standardize=TRUE,penaltyx=rv_cca_fit_perm$bestpenaltyx, K=4,
                 penaltyz=rv_cca_fit_perm$bestpenaltyz,v=rv_cca_fit_perm$v.init)

rv_radiomics_u1 = data.frame("feature"= colnames(x_fil)[which(rv_cca_fit$u[,1] != 0)],              #1st canonical vector of Radiomics Features
                             "coef" = as.numeric(rv_cca_fit$u[,1][which(rv_cca_fit$u[,1] != 0)]))   
rv_viral_v1 = data.frame("feature"= colnames(z_fil)[which(rv_cca_fit$v[,1] != 0)],                  #1st canonical vector of Viral Genome Word2Vec Encoding
                             "coef" = as.numeric(rv_cca_fit$v[,1][which(rv_cca_fit$v[,1] != 0)]))

rv_radiomics_u2 = data.frame("feature"= colnames(x_fil)[which(rv_cca_fit$u[,2] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$u[,2][which(rv_cca_fit$u[,2] != 0)]))
rv_viral_v2 = data.frame("feature"= colnames(z_fil)[which(rv_cca_fit$v[,2] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$v[,2][which(rv_cca_fit$v[,2] != 0)]))

rv_radiomics_u3 = data.frame("feature"= colnames(x_fil)[which(rv_cca_fit$u[,3] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$u[,3][which(rv_cca_fit$u[,3] != 0)]))
rv_viral_v3 = data.frame("feature"= colnames(z_fil)[which(rv_cca_fit$v[,3] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$v[,3][which(rv_cca_fit$v[,3] != 0)]))

rv_radiomics_u4 = data.frame("feature"= colnames(x_fil)[which(rv_cca_fit$u[,4] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$u[,4][which(rv_cca_fit$u[,4] != 0)]))
rv_viral_v4 = data.frame("feature"= colnames(z_fil)[which(rv_cca_fit$v[,4] != 0)], 
                             "coef" = as.numeric(rv_cca_fit$v[,4][which(rv_cca_fit$v[,4] != 0)]))
