import pandas as pd
import scdcdm
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as pl
import numpy as np
import arviz as az

from scdcdm.util import comp_ana as mod
from scdcdm.util import cell_composition_data as dat


#set data path
data_path = './data/'

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
tmp = meta.groupby(['Subject', 'Location', 'Replicate', 'Health','Epithelial_test'])['Subject'].count().unstack().fillna(0)
tmp = tmp.reset_index()
tmp_epi = tmp.loc[tmp['Location']=='Epi']

test_data = dat.from_pandas(tmp_epi, covariate_columns=["Health", 'Subject', 'Location', 'Replicate'])
print('Test data created')
#model set up (use three different references: 
for baseline in ['M cells', 'Enteroendocrine']: #Glia
    model_gut = mod.CompositionalAnalysis(test_data, formula="Health", baseline_index=baseline)
    print(baseline + ' cells.')
    for n_iter in np.array([20000, 40000, 80000, 150000], dtype='int32'):
        print(str(n_iter) + ' iterations.')
        sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
        sig_effects = sim_results.effect_df.loc[sim_results.effect_df["Final Parameter"] != 0]
        az.plot_density(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_' + baseline + '_' + str(n_iter) + '_density.pdf')
        az.plot_trace(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_' + baseline + '_' + str(n_iter) + '_trace.pdf')
        print(sig_effects)  
        #prepare result table to save to file
        sig_effects = sig_effects.reset_index()
        sig_effects['Health'] = sig_effects['Covariate'].astype('category')
        sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
        sig_effects['Epithelial_test'] = sig_effects['Cell Type']
        sig_effects.to_csv(data_path + 'results_' + baseline + '_' + str(n_iter) + '.csv')
        #sim_results.save(path_to_file = data_path + 'results_' + baseline + '_' + str(n_iter) +'.pkl')


#perform pairwise tests
print('Pairwise tests')

#Healthy vs Non-Inflamed
test_data_hn = test_data[test_data.obs['Health'] != 'Inflamed'].copy()
print('Test data created')
#model set up (use three different references: 
for baseline in ['Enterocyte Progenitors', 'TA 2', 'Cycling TA']: #Glia
    model_gut = mod.CompositionalAnalysis(test_data, formula="Health", baseline_index=baseline)
    print(baseline + ' cells.')
    for n_iter in np.array([20000, 40000, 80000, 150000], dtype='int32'):
        print(str(n_iter) + ' iterations.')
        sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
        sig_effects = sim_results.effect_df.loc[sim_results.effect_df["Final Parameter"] != 0]
        az.plot_density(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_hn_' + baseline + '_' + str(n_iter) + '_density.pdf')
        az.plot_trace(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_hn_' + baseline + '_' + str(n_iter) + '_trace.pdf')
        print(sig_effects)
        #prepare result table to save to file
        sig_effects = sig_effects.reset_index()
        sig_effects['Health'] = sig_effects['Covariate'].astype('category')
        sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
        sig_effects['Epithelial_test'] = sig_effects['Cell Type']
        sig_effects.to_csv(data_path + 'results_hn_' + baseline + '_' + str(n_iter) + '.csv')


#Healthy vs inflamed
test_data_hn = test_data[test_data.obs['Health'] != 'Non-inflamed'].copy()
print('Test data created')
#model set up (use three different references: 
for baseline in ['Goblet', 'Enteroendocrine']: #Glia
    model_gut = mod.CompositionalAnalysis(test_data, formula="Health", baseline_index=baseline)
    print(baseline + ' cells.')
    for n_iter in np.array([20000, 40000, 80000, 150000], dtype='int32'):
        print(str(n_iter) + ' iterations.')
        sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
        sig_effects = sim_results.effect_df.loc[sim_results.effect_df["Final Parameter"] != 0]
        az.plot_density(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_hi_' + baseline + '_' + str(n_iter) + '_density.pdf')
        az.plot_trace(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_hi_' + baseline + '_' + str(n_iter) + '_trace.pdf')
        print(sig_effects)
        #prepare result table to save to file
        sig_effects = sig_effects.reset_index()
        sig_effects['Health'] = sig_effects['Covariate'].astype('category')
        sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
        sig_effects['Epithelial_test'] = sig_effects['Cell Type']
        sig_effects.to_csv(data_path + 'results_hi_' + baseline + '_' + str(n_iter) + '.csv')


#non-inflamed vs inflamed
test_data_hn = test_data[test_data.obs['Health'] != 'Healthy'].copy()
print('Test data created')
#model set up (use three different references: 
for baseline in ['Tuft', 'M cells']: #Glia
    model_gut = mod.CompositionalAnalysis(test_data, formula="Health", baseline_index=baseline)
    print(baseline + ' cells.')
    for n_iter in np.array([20000, 40000, 80000, 150000], dtype='int32'):
        print(str(n_iter) + ' iterations.')
        sim_results = model_gut.sample_hmc(n_burnin=10000, num_results=n_iter)
        sig_effects = sim_results.effect_df.loc[sim_results.effect_df["Final Parameter"] != 0]
        az.plot_density(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_ni_' + baseline + '_' + str(n_iter) + '_density.pdf')
        az.plot_trace(sim_results, var_names="beta")
        pl.savefig(data_path + 'results_ni_' + baseline + '_' + str(n_iter) + '_trace.pdf')
        print(sig_effects)
        #prepare result table to save to file
        sig_effects = sig_effects.reset_index()
        sig_effects['Health'] = sig_effects['Covariate'].astype('category')
        sig_effects['Health'] = sig_effects['Health'].cat.rename_categories({'Health[T.Inflamed]': 'Inflamed',
                                                                           'Health[T.Non-inflamed]': 'Non-inflamed'})
        sig_effects['Epithelial_test'] = sig_effects['Cell Type']
        sig_effects.to_csv(data_path + 'results_ni_' + baseline + '_' + str(n_iter) + '.csv')
