#Date created: 12/10/2020
#Date modified: 19/11/2020
#Author created: Maren Büttner
#Author modified:
#Description: Differences in abundance test with Dirichlet Regression

#Project: scCODA paper

library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(DirichletReg)

# Load data
data_dir <- './data/'

metadata_table <- read_csv(paste0(data_dir, 'meta_processed.csv'))

# Create count table
#
# 
library(dplyr)

#prepare count table and proportions for Epi
df_epi <- metadata_table[metadata_table$Location=='Epi',]
df_epi$test <- df_epi$Cluster
df_epi$test[df_epi$Major_l1 %in% c('Fibroblasts', 'Immune')] = 'Non-Epithelial' 
df_epi <- df_epi %>% count(Subject, Replicate, Health, test)
df_epi <- df_epi %>% group_by(Subject, Health, Replicate) %>% 
  mutate(percent = n/sum(n))
#prepare count table and proportions for LP
df_lp <- metadata_table[metadata_table$Location=='LP',]
df_lp$test <- df_lp$Cluster
df_lp$test[df_lp$Major_l1 %in% c('Epithelial')] = 'Epithelial' 
df_lp <- df_lp %>% count(Subject, Replicate, Health, test)
df_lp <- df_lp %>% group_by(Subject, Health, Replicate) %>% 
  mutate(percent = n/sum(n))

#prepare pairwise tests 
pairs_list <- data.frame('state_1' = c("Healthy", "Healthy", "Non-inflamed"), 
                         'state_2' = c("Non-inflamed", "Inflamed", "Inflamed"))


#####################################
#Prepare data for test with DirichletReg
df_tmp <- df_lp #df_epi
location <- 'LP' #'Epi'
# Calculate regression
for (idx in 1:dim(pairs_list)[1]){
  state_1 <- as.character(pairs_list$state_1[idx])
  state_2 <- as.character(pairs_list$state_2[idx])
  print(paste0('Test ', state_1, ' vs ', state_2))
  
  select_ID <-df_tmp[['Health']] %in% c(state_1, state_2)
  metaDF_tmp <-  df_tmp[select_ID,]
  
  rowname_id <- ''
  for (idx in 1:dim(metaDF_tmp)[1]){
    rowname_id[idx] <- paste0(metaDF_tmp$Subject[idx], '_', metaDF_tmp$Replicate[idx], '_' , metaDF_tmp$Health[idx])
  }
  metaDF_tmp$rowname_ID <- rowname_id
  
  unique_names <- strsplit(unique(rowname_id), '_')
  tmp_cov <- data.frame('Subject'= sapply(unique_names, '[', 1), 
             'Replicate' = sapply(unique_names, '[', 2), 
             'Health' = sapply(unique_names, '[', 3))
  #prepare data
  #reshape data table into matrix
  dummy_mat <- matrix(0, nrow=length(unique(rowname_id)), ncol = length(unique(metaDF_tmp$test)))
  rownames(dummy_mat) = unique(rowname_id)
  colnames(dummy_mat) = unique(metaDF_tmp$test)
  
  for (row_id in 1:dim(dummy_mat)[1]){
    for (col_id in 1:dim(dummy_mat)[2]){
      n_tmp <- metaDF_tmp$percent[metaDF_tmp$rowname_ID == rownames(dummy_mat)[row_id] & 
                              metaDF_tmp$test == colnames(dummy_mat)[col_id]  
                              ]
      if (!length(n_tmp)==0){
        dummy_mat[row_id, col_id] <- n_tmp
      }
      }
  } 
  
  test <- as.data.frame(dummy_mat)
  test$counts <- DR_data(test)
  data = cbind(test, tmp_cov)
  
  #run fit
  fit = DirichReg(counts ~ Health, data)
  
  # Get p-values
  u = summary(fit)
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('Health', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  fit$pvals = pvals
  
  write.csv(x = fit$pvals, file = paste0(data_dir, 'DirichletReg_DA_test_', location, '_', state_1 ,'_vs_', state_2 ,'.csv'))
  
}


##########################
# for Dirichlet Regression, correct all p-values
#Epi
######
all_files = list.files(path=data_dir, pattern='DirichletReg_DA_test_Epi_')
all_pvals_coarse = data.frame()
for (files in all_files){
    tmp <- read.csv(paste0(data_dir, files), row.names = 'X')
    row.names(tmp) <- gsub('_', ' ', 
                           gsub('.csv', '', 
                                gsub('DirichletReg_DA_test_Epi_', '', files)))
    all_pvals_coarse <- rbind(all_pvals_coarse, tmp)
 
} 
#adjust pvalues (turn them into a vector by rows first)

pvals_coarse_adjust <-p.adjust(as.vector(t(all_pvals_coarse)), method='BH')
pvals_coarse <- data.frame(matrix(as.vector(t(pvals_coarse_adjust)), nrow = 3, byrow = TRUE), 
                           row.names = row.names(all_pvals_coarse)
)
colnames(pvals_coarse) <- colnames(all_pvals_coarse)

write.csv(x = pvals_coarse, file = paste0(data_dir , 'DirichletReg_DA_test_BH_adjusted_pvals_Epi.csv'))

# for Dirichlet Regression, correct all p-values
# LP
#######
all_files = list.files(path=data_dir, pattern='DirichletReg_DA_test_LP_')
all_pvals_coarse = data.frame()
for (files in all_files){
  tmp <- read.csv(paste0(data_dir, files), row.names = 'X')
  row.names(tmp) <- gsub('_', ' ', 
                         gsub('.csv', '', 
                              gsub('DirichletReg_DA_test_Epi_', '', files)))
  all_pvals_coarse <- rbind(all_pvals_coarse, tmp)
  
} 
#adjust pvalues (turn them into a vector by rows first)

pvals_coarse_adjust <-p.adjust(as.vector(t(all_pvals_coarse)), method='BH')
pvals_coarse <- data.frame(matrix(as.vector(t(pvals_coarse_adjust)), nrow = 3, byrow = TRUE), 
                           row.names = row.names(all_pvals_coarse)
)
colnames(pvals_coarse) <- colnames(all_pvals_coarse)

write.csv(x = pvals_coarse, file = paste0(data_dir , 'DirichletReg_DA_test_BH_adjusted_pvals_LP.csv'))
