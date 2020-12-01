import pandas as pd
import scdcdm
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as pl
import numpy as np
import arviz as az
import random as rn

from scdcdm.util import comp_ana as mod
from scdcdm.util import cell_composition_data as dat

#set data path
data_path = './data/'
result_path = data_path + 'subsample_donor_list/'
result_table_path = data_path + 'results_subsample/'

#read data
meta = pd.read_csv(data_path + 'meta_processed.csv')
print('Data read')
#add test category (with Epithelial cells as rest class
meta['Cluster'] = meta['Cluster'].astype('category')
meta['Epithelial_test'] = meta['Cluster'].cat.add_categories(['Non-Epithelial'])
meta['Epithelial_test'][np.in1d(meta['Major_l1'], ['Fibroblasts', 'Immune'])] = 'Non-Epithelial'
meta['Epithelial_test'] = meta['Epithelial_test'].cat.remove_unused_categories()
meta['Epithelial_test'] = meta['Epithelial_test'].astype(str)

#process data to select tests
#sum up the replicates as well
tmp = meta.groupby(['Subject', 'Location','Replicate', 'Health','Epithelial_test'])['Subject'].count().unstack().fillna(0)
tmp = tmp.reset_index()
tmp_epi = tmp.loc[tmp['Location']=='Epi']

donor_list = np.unique(tmp['Subject'])
n_iter = 150000

#sample the number of items to remove from list
#Subsample, remove one from the Healthy and the UC patients, respectively.
donor_healthy = np.unique(tmp.loc[tmp['Health']=='Healthy']['Subject'])
donor_UC = np.unique(tmp.loc[tmp['Health']!='Healthy']['Subject'])

for k in range(6,11,1):
    with open(result_path + "Epi_subsample_" + str(2*k) + ".txt", "a") as myfile:
        for idx in range(5):
            remove_donor = np.concatenate([rn.sample(set(donor_healthy), k=k), rn.sample(set(donor_UC), k=k)])
            #prepare for saving
            tmp_name = str(remove_donor).split('[')[1].split(']')[0]
            #Save for the use in R.
            myfile.write(tmp_name + '\n')
            meta_tmp = tmp_epi.loc[np.invert(np.in1d(tmp_epi['Subject'], remove_donor))]

            test_data = dat.from_pandas(meta_tmp, covariate_columns=["Health", 'Subject','Replicate', 'Location'])
            print('Add pseudocount')
            test_data.X = test_data.X +0.5
            print('Test data created')
            #model set up:
            baseline = 'Non-Epithelial'
            


            #perform pairwise tests
            print('Pairwise tests')
            #Healthy vs Non-Inflamed
            test_data_hn = test_data[test_data.obs['Health'] != 'Inflamed'].copy()
            print('Test data created')
            #model set up (use three different references: 
            
            model_gut = mod.CompositionalAnalysis(test_data_hn, formula="Health", baseline_index=baseline)

            print(str(n_iter) + ' iterations.')
            sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
            sig_effects = sim_results.effect_df
            az.plot_density(sim_results, var_names="beta")
            pl.savefig(result_table_path + 'results_pseudocount_subsample_hn_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_density.pdf')
            az.plot_trace(sim_results, var_names="beta")
            pl.savefig(result_table_path + 'results_pseudocount_subsample_hn_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_trace.pdf')
            print(sig_effects)
            #prepare result table to save to file
            sig_effects = sig_effects.reset_index()
            sig_effects['Health'] = sig_effects['Covariate'].astype('category')
            sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
            sig_effects['Epithelial_test'] = sig_effects['Cell Type']
            sig_effects.to_csv(result_table_path + 'results_pseudocount_subsample_hn_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '.csv')
           
            #Healthy vs inflamed
            test_data_hn = test_data[test_data.obs['Health'] != 'Non-inflamed'].copy()
            print('Test data created')
            #model set up (use three different references: 

            model_gut = mod.CompositionalAnalysis(test_data_hn, formula="Health", baseline_index=baseline)

            print(str(n_iter) + ' iterations.')
            sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
            sig_effects = sim_results.effect_df
            az.plot_density(sim_results, var_names="beta")
            pl.savefig(result_table_path + 'results_pseudocount_subsample_hi_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_density.pdf')
            az.plot_trace(sim_results, var_names="beta")
            pl.savefig(result_table_path + 'results_pseudocount_subsample_hi_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_trace.pdf')
            print(sig_effects)
            #prepare result table to save to file
            sig_effects = sig_effects.reset_index()
            sig_effects['Health'] = sig_effects['Covariate'].astype('category')
            sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
            sig_effects['Epithelial_test'] = sig_effects['Cell Type']
            sig_effects.to_csv(result_table_path + 'results_pseudocount_subsample_hi_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '.csv')
                
             #non-inflamed vs inflamed
            test_data_hn = test_data[test_data.obs['Health'] != 'Healthy'].copy()
            print('Test data created')
            #model set up: 

            model_gut = mod.CompositionalAnalysis(test_data_hn, formula="Health", baseline_index=baseline)
            print(baseline + ' cells.')
            print(str(n_iter) + ' iterations.')
            sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
            sig_effects = sim_results.effect_df
            az.plot_density(sim_results, var_names="beta")
            pl.savefig(result_table_path+ 'results_pseudocount_subsample_ni_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_density.pdf')
            az.plot_trace(sim_results, var_names="beta")
            pl.savefig(result_table_path + 'results_pseudocount_subsample_ni_' + baseline + '_' + str(n_iter) +'_' + str(2*k) + '_' + str(idx) + '_trace.pdf')
            print(sig_effects)
            #prepare result table to save to file
            sig_effects = sig_effects.reset_index()
            sig_effects['Health'] = sig_effects['Covariate'].astype('category')
            sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                       'Health[T.Non-inflamed]': 'Non-inflamed'})
            sig_effects['Epithelial_test'] = sig_effects['Cell Type']
            sig_effects.to_csv(result_table_path + 'results_pseudocount_subsample_ni_' + baseline + '_' + str(n_iter) + '_' + str(2*k) + '_' + str(idx) +'.csv')
        
