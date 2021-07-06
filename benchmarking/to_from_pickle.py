# This script contains utilities to transfer the pickle files from the benchmark studies to universally readable files (csv, h5ad)#

import pandas as pd
import pickle as pkl
import os
import anndata as ad
import warnings
import re
import numpy as np

warnings.filterwarnings("ignore", category=FutureWarning)


def benchmark_datasets_to_csv(pickle_path, uni_path):
    """
    Converts all .pkl files from generate_data to a csv and a h5 file in a new directory for long-term storage

    :param pickle_path: path to directory where pickled files (output of generate_data) are located
    :param uni_path: path to directory where data in csv and h5 format should be located
    :return:
    """

    k = 0

    data_names = os.listdir(pickle_path)

    if len(data_names) > 0:

        all_parameters = []
        all_datasets = {}

        for name in data_names:
            with open(pickle_path + name, "rb") as f:
                data = pkl.load(f)

            all_parameters.append(data["parameters"])

            for d in data["datasets"]:
                all_datasets[k] = d
                k += 1

        concat = ad.concat(all_datasets, label="dataset_no")
        params = pd.concat(all_parameters, ignore_index=True)

        params.to_csv(uni_path + "generation_parameters")
        concat.write_h5ad(uni_path + "generated_data")

    else:
        print(f"no files in {pickle_path}. skipping...")



def benchmark_datasets_to_pickle(uni_path, pickle_path):
    """
    Converts a csv and a h5 file to .pkl files that are identical to the ones from from generate_data.
    Reverse operation to benchmark_datasets_to_csv

    :param uni_path: path to directory where data in csv and h5 format are located
    :param pickle_path: path to directory where pickled files (equal to output of generate_data) should be located
    :return:
    """

    parameters = pd.read_csv(uni_path + "generation_parameters", index_col=0)

    adata = ad.read_h5ad(uni_path + "generated_data")
    datasets = []
    for i in adata.obs["dataset_no"].unique():
        datasets.append(adata[adata.obs["dataset_no"] == i].copy())

    batch_size = parameters.groupby(["Base", "Increase", "n_controls", "n_cases"]).size().iloc[0]

    i_s = []
    k = 0
    for i in range(parameters.shape[0]):

        i_s.append(i)

        if (i+1) % batch_size == 0:

            print(i)
            out = {"parameters": parameters.iloc[i_s,:].reset_index(drop=True),
                   "datasets": [datasets[x] for x in i_s]}

            with open(pickle_path + "data_" + str(k), "wb") as f:
                pkl.dump(out, f)

            k += 1
            i_s = []


def benchmark_results_to_csv(pickle_path, uni_path):
    """
    Converts all .pkl files from a benchmark to a csv and h5 files in a new directory for long-term storage

    :param pickle_path: path to directory where pickled files (output of benchmark) are located
    :param uni_path: path to directory where data in csv and h5 format should be located
    :return:
    """

    data_names = os.listdir(pickle_path)

    effects = {}

    if len(data_names) > 0:

        all_results = []
        sccoda_effects = []

        for name in data_names:
            if not name == ".":
                with open(pickle_path + name, "rb") as f:
                    temp = pkl.load(f)


                if name.startswith("scCODA"):
                    all_results.append(temp["results"])
                    sccoda_effects = sccoda_effects + temp["effects"]

                elif isinstance(temp, dict):
                    all_results.append(temp["results"])

                else:
                    all_results.append(temp)


        all_results = pd.concat(all_results)
        all_results.to_csv(uni_path + "benchmark_results")

        sccoda_conc = pd.concat(sccoda_effects, keys=list(range(len(sccoda_effects))), names=["dataset", "Covariate", "Cell Type"])
        sccoda_conc.to_csv(uni_path + "sccoda_effects")

    else:
        print(f"no files in {pickle_path}. skipping...")


def model_comparison_to_csv(pickle_path, uni_path):

    test_names = pd.unique([re.sub(r"_results[_]*[0-9]*.pkl", "", x) for x in os.listdir(pickle_path)])
    test_names = np.setdiff1d(test_names, ".DS_Store")

    effects = {}
    results = []
    for name in test_names:
        eff = []
        for file in os.listdir(pickle_path):
            if file.startswith(name):
                with open(pickle_path + file, "rb") as f:
                    temp = pkl.load(f)
                if name == "ancom":
                    for x in temp["effects"]:
                        if type(x[0]) == bool:
                            t = pd.DataFrame({
                                "W": np.zeros(len(x)),
                                "Reject null hypothesis": x
                            })
                            eff.append(t)
                        else:
                            eff.append(x[0])
                elif name == "scdc":
                    for x in temp["effects"]:
                        eff.append(x.iloc[int(len(x) / 2):, :])
                else:
                    for x in temp["effects"]:
                        eff.append(x)
                results.append(temp["results"])

        effects[name] = pd.concat(eff)
        pd.concat(eff).to_csv(uni_path + f"{name}_effects")

    results_ = pd.concat(results)
    results_.to_csv(uni_path + "benchmark_results")
