#Date created: 24/11/2020
#Date modified:  30/11/2020
#Author created: Maren Büttner
#Author modified:
#Description: Differences in abundance test with Dirichlet Regression upon subsampling

#Project: scCODA paper

library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(DirichletReg)

#get date
today = format(Sys.Date(), '%y%m%d')

# Load data
data_dir <- './data/'

result_dir <-'./results/subsampling/'
result_dir2 <-'./results/'

metadata_table <- read_csv(paste0(data_dir, 'meta_processed.csv'))

#Load subsampling info:
subsample_dir <- './results/subsample_donor_list/'
#get subsampling design
subsamples = seq(12, 20, 2)
iterations = 1:5

#test 
subsample = 12
subsample_ref <- read.csv(paste0(subsample_dir, 'Epi_subsample_', subsample, '.txt'), 
                          header = FALSE, quote = "\'", stringsAsFactors = FALSE)


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


#prepare pairwise tests 
pairs_list <- data.frame('state_1' = c("Healthy", "Healthy", "Non-inflamed"), 
                         'state_2' = c("Non-inflamed", "Inflamed", "Inflamed"))


#####################################
#Prepare data for test with DirichletReg
df_tmp <- df_epi
location <- 'Epi'

for (subsample in subsamples){
  #read donor list to exclude
  subsample_ref <- read.csv(paste0(subsample_dir, 'Epi_subsample_', subsample, '.txt'), 
                            header = FALSE, quote = "\'", stringsAsFactors = FALSE)
   
  dim_ref = dim(subsample_ref)
  #prepare subsample_ref, because some files had line breaks
  if (dim_ref[2]!=subsample){
    list_ref <- list()
    len_list <- list()
    for(idx in 1:dim_ref[1]) {
      #split the read table into the respective samples as array of strings
      subsample_ref[idx,] <- trimws(subsample_ref[idx,], which = "left")
      list_ref[[idx]] <- strsplit(subsample_ref[idx,], ' ')[[1]]
      len_list[idx] <- length(list_ref[[idx]])
    }
    subsample_ref_new <- matrix(ncol=subsample, nrow=0)
    for (idx in 1:dim_ref[1]){
      if(len_list[idx]==subsample){
        subsample_ref_new <- rbind(subsample_ref_new, list_ref[[idx]])
      }else{ #concatenate two lines
        if (idx ==dim_ref[1]){ next } #skip last line
        tmp <- c(list_ref[[idx]], list_ref[[idx+1]])
        if (length(tmp)==subsample){
          subsample_ref_new <- rbind(subsample_ref_new, tmp)
        }
      }
    }
  }
  subsample_ref <- subsample_ref_new
  
  for (k in iterations){
    
    #trim trailing whitespace
    ref_tmp <- trimws(subsample_ref[k,], which = "left")
    #subsample input data
    df_sub_tmp <- df_tmp %>% filter(!Subject %in% ref_tmp)

    
    # Calculate regression
    for (idx in 1:dim(pairs_list)[1]){
      state_1 <- as.character(pairs_list$state_1[idx])
      state_2 <- as.character(pairs_list$state_2[idx])
      print(paste0('Test ', state_1, ' vs ', state_2))
      
      select_ID <-df_sub_tmp[['Health']] %in% c(state_1, state_2)
      metaDF_tmp <-  df_sub_tmp[select_ID,]
      
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
      
      write.csv(x = fit$pvals, file = paste0(result_dir, 'DirichletReg_DA_test_', 
                                             location, '_subsample_', 
                                             subsample, '_', k, '_',
                                             state_1 ,'_vs_', state_2 ,'.csv'))
      
    }
    
    
    # for Dirichlet Regression, correct all p-values
    prefix =paste0('DirichletReg_DA_test_Epi_subsample_', subsample, '_', k, '_')
    all_files = list.files(path=result_dir, 
                           pattern=prefix)
    all_pvals_coarse = data.frame()
    for (files in all_files){
      tmp <- read.csv(paste0(result_dir, files), row.names = 'X')
      row.names(tmp) <- gsub('_', ' ', 
                             gsub('.csv', '', 
                                  gsub(prefix, '', files)))
      all_pvals_coarse <- rbind(all_pvals_coarse, tmp)
      
    } 
    #adjust pvalues (turn them into a vector by rows first)
    
    pvals_coarse_adjust <-p.adjust(as.vector(t(all_pvals_coarse)), method='BH')
    pvals_coarse <- data.frame(matrix(as.vector(t(pvals_coarse_adjust)), nrow = 3, byrow = TRUE), 
                               row.names = row.names(all_pvals_coarse)
    )
    colnames(pvals_coarse) <- colnames(all_pvals_coarse)
    
    write.csv(x = pvals_coarse, file = paste0(data_dir ,
                                      'DirichletReg_DA_test_BH_adjusted_pvals_Epi_subsample_',
                                      subsample, '_', k, 
                                      '.csv'))
    
  }
}


##################
# compute MCC for all results
# use significance values from the full dataset as reference
#
library(mltools)

#read reference data
ref_DA <- read.csv(paste0(data_dir, 'DirichletReg_DA_test_BH_adjusted_pvals_Epi.csv')) 
#set row names and remove first column
row.names(ref_DA) <- ref_DA$X
ref_DA <- ref_DA[,-1]

#create significance indicator
ref_signif <- ref_DA<0.05

#create dataframe for MCC results
results_mcc <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("MCC", "Tested levels", "Subsample", "Iteration")
colnames(results_mcc) <- x

#run comparison with reference
for (subsample in subsamples){
  for (k in iterations){
    test_data <- read.csv(paste0(data_dir ,
                                  'DirichletReg_DA_test_BH_adjusted_pvals_Epi_subsample_',
                                  subsample, '_', k, 
                                  '.csv'))
    #set row names and remove first column
    row.names(test_data) <- test_data$X
    test_data <- test_data[,-1]
    #create significance indicator
    test_signif <- test_data < 0.05
    for (level in 1:3){
      tmp_mcc <-mcc (preds = test_signif[level, ], actuals = ref_signif[level,])
      results_mcc[nrow(results_mcc) + 1,] = c(tmp_mcc, 
                                              rownames(ref_signif)[level], 
                                              subsample,
                                              k
                                              )
      
    }
  }
}

#save result MCC to file.
write.csv(x = results_mcc, file = paste0(result_dir2, today, 'DirichletReg_DA_mcc_subsampling_epi.csv'))

