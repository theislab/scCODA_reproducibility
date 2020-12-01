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
meta['LP_test'] = meta['Cluster'].cat.add_categories(['Epithelial'])
meta['LP_test'][np.in1d(meta['Major_l1'], ['Epithelial'])] = 'Epithelial'
meta['LP_test'] = meta['LP_test'].cat.remove_unused_categories()
meta['LP_test'] = meta['LP_test'].astype(str)
#process data to select tests
tmp = meta.groupby(['Subject', 'Location', 'Replicate', 'Health','LP_test'])['Subject'].count().unstack().fillna(0)
tmp = tmp.reset_index()
tmp_lp = tmp.loc[tmp['Location']=='LP']

test_data = dat.from_pandas(tmp_lp, covariate_columns=["Health", 'Subject', 'Location', 'Replicate'])
print('Test data created')
#model set up (use three different references: Glia, Endothelial, and Cycling Monocytes
for baseline in ['Post-capillary Venules', 'Inflammatory Monocytes']: #Glia
    test_data_hn = test_data[test_data.obs['Health'] != 'Inflamed'].copy()
    model_gut = mod.CompositionalAnalysis(test_data_hn, formula="Health", baseline_index=baseline)
    print(baseline + ' cells.')
    for n_iter in np.array([200000, 400000, 800000], dtype='int32'):
        print(str(n_iter) + ' iterations.')
        sim_results = model_gut.sample_hmc(n_burnin=20000, num_results=n_iter)
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
        sig_effects['LP_test'] = sig_effects['Cell Type']
        sig_effects.to_csv(data_path + 'results_hn_' + baseline + '_' + str(n_iter) + '.csv')
        #sim_results.save(path_to_file = data_path + 'results_' + baseline + '_' + str(n_iter) +'.pkl')



