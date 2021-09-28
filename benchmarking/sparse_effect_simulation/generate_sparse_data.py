# Generates compositional datasets for benchmarks.

import numpy as np
import pickle as pkl
import pandas as pd
import itertools
import sys
import os

# Insert path to scCODA package for running on server
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/scCODA/')

from sccoda.util import data_generation as gen


def generate_compositional_datasets_sparse(n_cell_types, n_cells, n_samples,
                                           fct_base, fct_change,
                                           n_repetitions, sparsity_level, mode="absolute",
                                           write_path=None, file_name=None):
    """
    Generate compositional case-control data for all combinations of n_cell_types, n_cells, n_samples,
    fct_base, fct_change and save them to disk or return them.
    Datasets are always modeled such that there is an effect on the first cell type, while all other types are unaffected
    For each parameter combination, n_repetitions datasets are generated.

    Parameters
    ----------
    n_cell_types: list
        Number of cell types
    n_cells: list
        total number of cells per sample
    n_samples: list
        Number of samples per group. Each list enement must be of the type [n_controls, n_cases]
    fct_base: list
        Mean abundance for first cell type
    fct_change: list
        Change in first cell type between groups. See mode for details
    n_repetitions: int
        Number of repeated data generations for each parameter combination
    mode: str
        If "absolute", fct_change is interpreted as an absolute change in cell counts
        If "relative", fct_change is interpreted as a change relative to fct_base.
        fct_change(absolute) = fct_change(relative)*fct_base
    write_path: str
        Path to folder where files are written. If None, data is returned instead of written to disk
    file_name
        prefix for files. If None, data is returned instead of written to disk

    Returns
    -------
    If writing to disk is chosen, writes one pickled file per combination of parameters to disk.
    Each file contains n_repetitions datasets.
    They are structured a dict with "parameters" being a DataFrame that contains the generation parameters
    and "datasets" being a list of scCODA datasets (see the scCODA documentation for details)

    Otherwise, all generated datasets are returned in one dict with the same structure as described above
    """

    # Cenerate list of all parameter combinations

    if mode == "absolute":
        if n_cells == "adaptive":
            temp_params = list(itertools.product(n_cell_types, n_samples, sparsity_level, fct_base, fct_change))
            simulation_params = []
            for p in temp_params:
                p = list(p)
                n_cells = p[0] * p[2]
                p = [p[0], n_cells, p[1], p[2], p[3], p[4]]
                simulation_params.append(tuple(p))
        else:
            simulation_params = list(itertools.product(n_cell_types, n_cells, n_samples, sparsity_level, fct_base, fct_change))

    elif mode == "relative":
        # For relative mode, we need to calculate the change for each parameter set separately
        if n_cells == "adaptive":
            temp_params_ = list(itertools.product(n_cell_types, n_samples, sparsity_level, fct_base))
            temp_params = []
            for p in temp_params_:
                p = list(p)
                n_cells = p[0] * p[3]
                p = [p[0], n_cells, p[1], p[2], p[3]]
                temp_params.append(tuple(p))
        else:
            temp_params = list(itertools.product(n_cell_types, n_cells, n_samples, sparsity_level, fct_base))

        simulation_params = []
        for p in temp_params:
            p = list(p)

            # Balanced base generation
            if p[4] == "balanced":
                p[4] = p[1] * (1 / p[0])
            for c in fct_change:
                base = np.round(gen.counts_from_first(p[4], p[1], p[0]), 3)
                n_da = len(c)
                change = base[0] * np.array(c)

                # bugfix if da cell types make up all the cells
                if np.sum(base[:n_da] + change) >= p[1]:
                    change = change - (np.sum(base[:n_da] + change) - p[1]) - 1
                p_ = p.copy()
                p_.append(tuple(change))
                p_ = tuple(p_)
                simulation_params.append(p_)
                print(p_)

    else:
        raise ValueError("Wrong mode specified!")

    # Initialize output components
    out_parameters = pd.DataFrame(columns=['n_cell_types', 'n_cells',
                                           'n_controls', 'n_cases', 'sparsity_level',
                                           'Base', 'Increase',
                                           'b_true', 'w_true'])
    out_datasets = []

    i = 1

    # iterate over all combinations of parameters
    for n_cell_types_, n_cells_, n_samples_, sparsity_level_, fct_base_, fct_change_ in simulation_params:

        # initialize parameter df
        parameters = pd.DataFrame(columns=['n_cell_types', 'n_cells',
                                           'n_controls', 'n_cases', 'sparsity_level',
                                           'Base', 'Increase', 'log-fold increase',
                                           'b_true', 'w_true'])
        datasets = []

        print(f"parameters {i}/{len(simulation_params)}")
        if type(fct_change_) == tuple:
            change_ = np.pad(fct_change_, (0, n_cell_types_ - len(fct_change_)))
        else:
            change_ = np.pad(fct_change_, (0, n_cell_types_ - 1))

        # calculate b for use in gen.generate_case_control, normalize
        b = np.round(gen.counts_from_first(fct_base_, n_cells_, n_cell_types_), 3)
        b_t = np.round(np.log(b / n_cells_), 3)
        # calculate w for use in gen.generate_case_control
        _, w = np.round(gen.b_w_from_abs_change(b, change_, n_cells_), 3)

        # Generate n_repetitions datasets, add to parameter df and dataset list
        for j in range(n_repetitions):
            sigma = np.identity(n_cell_types_) * 0.05
            n_unaff = int((1-sparsity_level_)*n_samples_[1])
            temp_data = gen.generate_case_control(cases=1, K=n_cell_types_,
                                                  n_total=n_cells_, n_samples=[n_samples_[0] + n_unaff, n_samples_[1] - n_unaff],
                                                  b_true=b_t, w_true=[w],
                                                  sigma=sigma)
            temp_data.obs["x_0"] = np.repeat([0,1], n_samples_)
            datasets.append(temp_data)

            params = [n_cell_types_, n_cells_,
                      n_samples_[0], n_samples_[1], sparsity_level_,
                      fct_base_, fct_change_,
                      b_t, [w]]
            parameters = parameters.append(dict(zip(parameters.columns, params)), ignore_index=True)

        # If writing to disk, do this now
        if write_path is not None and file_name is not None:
            write = {"datasets": datasets, "parameters": parameters}

            with open(write_path + "/" + file_name + "_" + str(i), "wb") as f:
                pkl.dump(write, f)
        # Else, add generated data and parameters to output objects
        else:
            out_datasets = out_datasets + datasets
            out_parameters = out_parameters.append(parameters)

        i += 1

    out = {"datasets": out_datasets, "parameters": out_parameters}

    return out

__name__="main"
if __name__ == "main":
    np.random.seed(1234)

    n_cell_types = [5]
    n_cells = [5000]
    n_samples = [[20, 20]]
    fct_base = [1, 100, 1000]
    fct_change = [500, 1000, 2000]
    sparsity_level = list(np.round(np.arange(0.05, 1.01, 0.05), 2))
    n_repetitions = 20

    write_path = os.path.realpath(
        "../other_data/sparsity_test/generated_datasets_model_sparsity_test/")
    file_name = "sparsity_data"

    comp_data = generate_compositional_datasets_sparse(n_cell_types=n_cell_types, n_cells=n_cells,
                                                          n_samples=n_samples, fct_base=fct_base, fct_change=fct_change,
                                                          sparsity_level=sparsity_level,
                                                          n_repetitions=n_repetitions, mode="absolute",
                                                       write_path=write_path, file_name=file_name)
