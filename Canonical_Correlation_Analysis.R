# Canonical correlation analysis of two data modalities.

library(PMA)
library(ggplot2)
library(dplyr)
library("data.table")
library(tidyr)

sparse_cca <- function(x,z){  # Assuming we have two data matrices, x and z, on the same set of n observations
  list_result <- list()
  x <- apply(as.matrix(x),2,as.numeric)
  z <- apply(as.matrix(z),2,as.numeric)
  X_ind = which(apply(x,2,var) > 0.001) 
  Z_ind = which(apply(z,2,var) > 0.001) 
  x_fil = x[,X_ind]
  z_fil = z[,Z_ind]
  
  set.seed(3)
  cca_fit_perm = CCA.permute(x_fil,z_fil,   # According the the number of canonical vectors desired, k parameter of the CCA function and the returned list can be adjusted.
                                typex="standard",
                                typez="standard",
                                nperms=25,trace=FALSE)
  
  cca_fit = CCA(x_fil,z_fil,typex = 'standard',typez='standard',standardize=TRUE,penaltyx=cca_fit_perm$bestpenaltyx,
                   penaltyz=cca_fit_perm$bestpenaltyz,v=cca_fit_perm$v.init)
  
  u1 = data.frame("feature"= colnames(x_fil)[which(cca_fit$u[,1] != 0)], 
                               "coef" = as.numeric(cca_fit$u[,1][which(cca_fit$u[,1] != 0)]))
  v1 = data.frame("feature"= colnames(z_fil)[which(cca_fit$v[,1] != 0)], 
                           "coef" = as.numeric(cca_fit$v[,1][which(cca_fit$v[,1] != 0)]))
  
  list_result[[length(list_result) + 1]] <- list(u1,v1)
  
  return(list_result) #First canonical vectors of relevant data modalities
  
}

