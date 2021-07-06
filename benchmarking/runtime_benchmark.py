import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
import seaborn as sns

from benchmarking.generate_data import generate_compositional_datasets
from sccoda.util import comp_ana as ana

sns.set_style("ticks")
sns.set_context("paper", font_scale=1.4)

#%%

np.random.seed(5678)

n_cell_types = np.arange(5, 55, 5).tolist()
n_cells = [100000]
n_samples = [[5, 5], [10, 10], [15, 15], [20, 20]]
fct_base = ["balanced"]
fct_change = [1]
n_repetitions = 20

comp_data = generate_compositional_datasets(n_cell_types=n_cell_types, n_cells=n_cells,
                                            n_samples=n_samples, fct_base=fct_base, fct_change=fct_change,
                                            n_repetitions=n_repetitions, mode="relative",
                                            write_path=None, file_name=None)

#%%
comp_results = []
n = 0

for test_data in comp_data["datasets"]:

    n += 1
    print(f"dataset {n}/{len(comp_data['datasets'])}")

    K = test_data.X.shape[1]
    test_model = ana.CompositionalAnalysis(test_data, "x_0", K-1)

    result_1 = test_model.sample_hmc(num_results=1000, num_burnin=500)
    result_2 = test_model.sample_hmc(num_results=2000, num_burnin=500)

    comp_results.append((result_1, result_2))


result_df = comp_data["parameters"].copy()

result_df["time_1000"] = [x[0].sampling_stats["duration"] for x in comp_results]
result_df["time_2000"] = [x[1].sampling_stats["duration"] for x in comp_results]

result_df["time_per_step"] = (result_df.loc[:, "time_2000"] - result_df.loc[:, "time_1000"]) / 1000

print(result_df)

result_save_path = "../zenodo_data/sccoda_benchmark_data/runtime_analysis/"
result_df.to_csv(result_save_path + "runtime_analysis_results_2")

#%%
result_save_path = os.path.abspath("../sccoda_benchmark_data/runtime_results/")

result_df = pd.read_csv(result_save_path + "/runtime_analysis_results_rev")
plot_df = result_df.copy()
plot_df.rename(columns={
    "n_cell_types": "Cell types",
    "time_per_step": "Time per step (s)",
    "n_controls": "Samples per group"
}, inplace=True)
print(plot_df)

palette = sns.color_palette("Blues", len(plot_df["Samples per group"].unique()))

fig, ax = plt.subplots(figsize=(12, 6))
sns.lineplot(
    data=plot_df,
    x="Cell types",
    y="Time per step (s)",
    hue="Samples per group",
    palette=palette,

)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title="Samples per group")
plt.tight_layout()

plt.savefig(result_save_path + "/runtime_lines_rev.png")
plt.savefig(result_save_path + "/runtime_lines_rev.svg")

plt.show()




#%%
pd.options.display.float_format = '{:10,.5f}'.format
pd.set_option('display.max_columns', 10)
pd.set_option('display.max_rows', None)

mean_df = result_df.groupby(["n_controls", "n_cell_types"]).agg(
    {"time_per_step": "mean"}
)
print(mean_df)

#%%
print(result_df.loc[
          (result_df["n_controls"].isin([5]))
          & (result_df["n_cell_types"].isin([50]))
      ])

#%%

print(np.max(result_df["time_per_step"]))
#%%

print(np.mean(result_df["time_per_step"]))
