# Utility functions for benchmarks on scCODA
import numpy as np
import pandas as pd
import pickle as pkl
import os
import sys
import patsy as pt
import rpy2

# For running on server
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/scCODA/')

from sccoda.util import data_generation as gen
# from sccoda.model import other_models as om
from sccoda.model import scCODA_model as model


def benchmark(data_path, save_path, models, benchmark_name="", server=False, keep_results=False, alpha=0.05):
    """
    Run a benchmark. All models in parameter "models" are applied to all datasets in "data_path".

    Parameters
    ----------
    data_path: str
        Path to folder where datasets are saved
    save_path: str
        Path to folder where results are saved
    models: list
        List of models to include in the benchmark
    benchmark_name: str
        prefix for meta files when executing on server.
    server: bool
        Execute on ICB server

    Returns
    -------
    Saves pickled DataFrames to disk that contain generation parameters, model name and  effect discovery results.

    """

    # Models that need batched execution due to calculation time (One data file at a time, also one results file per data file).
    # For all other models, only one results file with all results is generated
    batched_models = ["simple_dm", "scCODA", "scdc"]

    # Parameters for each model
    for model_name in models:
        print(model_name)

        if server:
            r_home = "/home/icb/johannes.ostner/anaconda3/lib/R"
            r_path = r"/home/icb/johannes.ostner/anaconda3/lib/R/bin"
        else:
            r_home = "/Library/Frameworks/R.framework/Resources"
            r_path = r"/Library/Frameworks/R.framework/Resources/bin"

        if model_name == "ALDEx2_alr":
            kwargs = {
                "r_home": r_home,
                "r_path": r_path,
                "method": "we.eBH",
                "mc_samples": 128,
                "denom": [5],
                "alpha": alpha,
                "fdr_correct": False,
                "keep_results": keep_results
            }

        elif model_name == "ALDEx2":
            kwargs = {
                "r_home": r_home,
                "r_path": r_path,
                "method": "we.eBH",
                "mc_samples": 128,
                "alpha": alpha,
                "fdr_correct": False,
                "keep_results": keep_results
            }

        elif model_name in ["simple_dm", "scCODA"]:
            kwargs = {"num_results": 20000,
                      "num_burnin": 5000,
                      "num_adapt_steps": 4000,
                      "keep_results": keep_results,
                      "alpha": alpha,
                      }

        elif model_name in ["alr_ttest", "alr_wilcoxon"]:
            kwargs = {"reference_cell_type": 4,
                      "alpha": alpha,
                      "fdr_correct": True,
                      "keep_results": keep_results
                      }
        elif model_name in ["Haber", "ttest", "clr_ttest"]:
            kwargs = {"alpha": alpha,
                      "fdr_correct": True,
                      "keep_results": keep_results
                      }
        elif model_name in ["scdc", "dirichreg", "BetaBinomial"]:
            kwargs = {
                "r_home": r_home,
                "r_path": r_path,
                "keep_results": keep_results,
                "alpha": alpha
            }
        elif model_name =="ancom":
            kwargs = {
                "keep_results": keep_results,
                "alpha": alpha
            }
        elif model_name == "ANCOMBC":
            kwargs = {
                "r_home": r_home,
                "r_path": r_path,
                "fdr_correct": False,
                "keep_results": keep_results,
                "alpha": alpha
            }
        elif model_name == "ANCOMBC_holm":
            kwargs = {
                "r_home": r_home,
                "r_path": r_path,
                "fdr_correct": False,
                "keep_results": keep_results,
                "method": "holm",
                "alpha": alpha
            }
        else:
            kwargs = {}

        # For each batched model, run one datafile at a time:
        if model_name in batched_models:
            num_files = len(os.listdir(data_path))

            for count in range(num_files):
                # On server, generate bash file and push to queue
                if server:
                    bash_loc = f"/home/icb/johannes.ostner/compositional_diff/benchmark_scripts/{benchmark_name}"
                    bash_name = f"{benchmark_name}_{model_name}_{count}"
                    script_location = "/home/icb/johannes.ostner/compositional_diff/benchmarking/model_one_job_batched.py"
                    arguments = [data_path, save_path, model_name, alpha, count, keep_results]

                    execute_on_server(bash_loc, bash_name, script_location, arguments)
                # Locally, just execute and save
                else:
                    file_name = os.listdir(data_path)[count]

                    if keep_results:
                        results, effects = model_on_one_datafile(data_path + file_name, model_name, **kwargs)
                        results = get_scores(results)
                        save = {"results": results, "effects": effects}

                    else:
                        results = model_on_one_datafile(data_path + file_name, model_name, **kwargs)

                        results = get_scores(results)
                        save = results

                    with open(save_path + model_name + "_results_" + str(count) + ".pkl", "wb") as f:
                        pkl.dump(save, f)
        # For unbatched models, run all datasets at once
        else:
            # On server, generate bash file and push to queue
            if server:
                bash_loc = f"/home/icb/johannes.ostner/compositional_diff/benchmark_scripts/{benchmark_name}"
                bash_name = f"{benchmark_name}_{model_name}"
                script_location = "/home/icb/johannes.ostner/compositional_diff/benchmarking/model_one_job.py"
                arguments = [data_path, save_path, model_name, alpha]

                execute_on_server(bash_loc, bash_name, script_location, arguments)
            # Locally, just execute and save
            else:
                results = model_all_datasets(data_path, model_name, **kwargs)

                if keep_results:
                    results["results"] = get_scores(results["results"])
                else:
                    results = get_scores(results)

                with open(save_path + model_name + "_results.pkl", "wb") as f:
                    pkl.dump(results, f)


def model_on_one_datafile(file_path, model_name, keep_results=False, *args, **kwargs):
    """
    Run a model on one datafile from generate_data.generate_compositional_datasets

    Parameters
    ----------
    file_path: str
        path to dataset
    model_name: str
        model name
    keep_results: bool
        Whether betas_df from the scCODA result should be returned as well (for threshold benchmark)
    args:
        Passed to the model execution function
    kwargs:
        Passed to the model execution function


    Returns
    -------
    fin_df: DataFrame
        contains generation parameters, model name and  effect discovery results

    """
    with open(file_path, "rb") as f:
        data = pkl.load(f)

    # Initialize output df
    fin_df = data["parameters"]
    fin_df["tp"] = 0
    fin_df["tn"] = 0
    fin_df["fp"] = 0
    fin_df["fn"] = 0
    fin_df["model"] = model_name

    # Parameters for model evaluation, not execution
    alpha = kwargs.pop("alpha", 0.05)
    tau = kwargs.pop("tau", 0.02)
    fdr_correct = kwargs.pop("fdr_correct", True)

    # initialize list for beta_dfs
    effect_dfs = []

    # Select right model to execute. For most models, just init model, run fit, run eval
    for d in range(len(data["datasets"])):

        model_completed = True
        K = fin_df["n_cell_types"].iloc[d]
        increase = fin_df["Increase"].iloc[d]
        try:
            inc = len(increase)
            ground_truth = np.array(increase)
            ground_truth = np.pad(ground_truth, (0, K - inc))
        except TypeError:
            ground_truth = np.zeros((K))
            ground_truth[0] = increase

        gt = (ground_truth != 0)

        if model_name in ["Haber"]:
            mod = om.HaberModel(data["datasets"][d])
            mod.fit_model()
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "ttest":
            mod = om.TTest(data["datasets"][d])
            mod.fit_model()
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "clr_ttest":
            mod = om.CLR_ttest(data["datasets"][d])
            mod.fit_model()
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name in ["ALDEx2", "ALDEx2_alr"]:
            mod = om.ALDEx2Model(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "alr_ttest":
            mod = om.ALRModel_ttest(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "alr_wilcoxon":
            mod = om.ALRModel_wilcoxon(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "ancom":
            mod = om.AncomModel(data["datasets"][d])
            mod.fit_model(alpha=alpha, tau=tau)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt)
            if keep_results:
                effect_dfs.append(mod.ancom_out)

        elif model_name in ["ANCOMBC", "ANCOMBC_holm"]:
            mod = om.ANCOMBCModel(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "BetaBinomial":
            mod = om.BetaBinomialModel(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "dirichreg":
            mod = om.DirichRegModel(data["datasets"][d])
            mod.fit_model(*args, **kwargs)
            tp, tn, fp, fn = mod.eval_model(ground_truth=gt, alpha=alpha, fdr_correct=fdr_correct)
            if keep_results:
                effect_dfs.append(pd.DataFrame({
                    "Cell Type": [str(x) for x in np.arange(K)],
                    "p value": mod.p_val,
                }))

        elif model_name == "simple_dm":
            # Build covariance matrix for simple dm and extract comp. data
            dat = data["datasets"][d]
            K = dat.X.shape[1]
            # Only one covariate
            formula = "x_0"

            cell_types = dat.var.index.to_list()

            # Get count data
            data_matrix = dat.X.astype("float32")

            # Build covariate matrix from R-like formula
            covariate_matrix = pt.dmatrix(formula, dat.obs)
            covariate_names = covariate_matrix.design_info.column_names[1:]
            covariate_matrix = covariate_matrix[:, 1:]

            # Init model. Baseline index is always the last cell type
            mod = om.SimpleModel(covariate_matrix=np.array(covariate_matrix), data_matrix=data_matrix,
                                 cell_types=cell_types, covariate_names=covariate_names, formula=formula,
                                 reference_cell_type=K-1)

            # Run HMC sampling
            result_temp = mod.sample_hmc(*args, **kwargs)
            alphas_df, betas_df = result_temp.summary_prepare(hdi_prob=1-alpha)

            # Effects are significant if 0 not in 95% HDI, set non-sig. effects to 0
            a_cols = [x for x in betas_df.columns if x.startswith("HDI")]
            final_betas = np.where((betas_df.loc[:, a_cols[0]] < 0) &
                                   (betas_df.loc[:, a_cols[1]] > 0),
                                   0,
                                   betas_df.loc[:, "Final Parameter"])

            # Compare found significances to ground truth (only first cell type significant)
            # Get true positives, true negatives, false positives, false negatives
            true_indices = np.where(gt == True)[0]
            false_indices = np.where(gt == False)[0]

            tp = sum(final_betas[true_indices] != 0)
            fn = sum(final_betas[true_indices] == 0)
            tn = sum(final_betas[false_indices] == 0)
            fp = sum(final_betas[false_indices] != 0)

            # append betas_df to output list
            if keep_results:
                effect_dfs.append(betas_df)

        elif model_name == "scCODA":
            # Build covariance matrix for scCODA and extract comp. data
            dat = data["datasets"][d]
            K = dat.X.shape[1]
            # Only one covariate
            formula = "x_0"

            cell_types = dat.var.index.to_list()

            # Get count data
            data_matrix = dat.X.astype("float32")

            # Build covariate matrix from R-like formula
            covariate_matrix = pt.dmatrix(formula, dat.obs)
            covariate_names = covariate_matrix.design_info.column_names[1:]
            covariate_matrix = covariate_matrix[:, 1:]

            # Init model. Baseline index is always the last cell type
            mod = model.scCODAModel(covariate_matrix=np.array(covariate_matrix), data_matrix=data_matrix,
                                      cell_types=cell_types, covariate_names=covariate_names, formula=formula,
                                      reference_cell_type=K-1)

            # Run HMC sampling, get results
            result_temp = mod.sample_hmc(*args, **kwargs)
            alphas_df, betas_df = result_temp.summary_prepare(est_fdr=alpha)
            final_betas = np.array(betas_df.loc[:, "Final Parameter"])

            # Compare found significances to ground truth (only first cell type significant)
            # Get true positives, true negatives, false positives, false negatives
            true_indices = np.where(gt == True)[0]
            false_indices = np.where(gt == False)[0]

            tp = sum(final_betas[true_indices] != 0)
            fn = sum(final_betas[true_indices] == 0)
            tn = sum(final_betas[false_indices] == 0)
            fp = sum(final_betas[false_indices] != 0)

            # append betas_df to output list
            if keep_results:
                effect_dfs.append(betas_df)

        elif model_name == "scdc":
            mod = om.scdney_model(data["datasets"][d])
            try:
                effect, (tp, tn, fp, fn) = mod.analyze(
                    r_home=kwargs["r_home"],
                    r_path=kwargs["r_path"],
                    ground_truth=gt,
                    alpha=alpha
                )
                if keep_results:
                    effect_dfs.append(effect)
            except:
                print("Error encountered!")
                tp, tn, fp, fn = (np.nan, np.nan, np.nan, np.nan)

        else:
            raise ValueError("Invalid model name specified!")

        fin_df.loc[d, "tp"] = tp
        fin_df.loc[d, "fp"] = fp
        fin_df.loc[d, "tn"] = tn
        fin_df.loc[d, "fn"] = fn

    if keep_results:
        return fin_df, effect_dfs
    else:
        return fin_df


def model_all_datasets(directory, model_name, keep_results=False, *args, **kwargs):
    """
    Runs one model on all dataset files in a directory

    Parameters
    ----------
    directory: str
        path to directory
    model_name: str
        model name
    args:
        Passed to the model execution function
    kwargs:
        Passed to the model execution function

    Returns
    -------
    results: DataFrame
        contains generation parameters, model name and  effect discovery results
    """

    file_names = os.listdir(directory)

    l = len(file_names)
    k = 1

    results = []

    # Run model on all files, coerce results
    if keep_results:
        effects = []
        for name in file_names:
            print(f"{k}/{l}")


            res, eff = model_on_one_datafile(directory + name, model_name, keep_results=True, *args, **kwargs)
            results.append(res)
            effects.extend(eff)
            k += 1

        results_df = pd.concat(results)
        out = {"results": results_df, "effects": effects}

        return out

    else:
        for name in file_names:
            print(f"{k}/{l}")

            res = model_on_one_datafile(directory+name, model_name, keep_results=False, *args, **kwargs)

            results.append(res)
            k += 1

        results_df = pd.concat(results)

        return results_df


def get_scores(df):
    """
    Calculates extended binary classification summary statistics, such as TPR, TNR, youden index, f1-score, MCC

    Parameters
    ----------
    df: DataFrame
        Must contain columns tp, tn, fp, fn

    Returns
    -------
    df: DataFrame
        Same df with added columns tpr, tnr, precision, accuracy, youden, f1_score, mcc
    """
    tp = df["tp"].astype("float64")
    tn = df["tn"].astype("float64")
    fp = df["fp"].astype("float64")
    fn = df["fn"].astype("float64")

    tpr = (tp / (tp + fn)).fillna(0)
    df["tpr"] = tpr
    tnr = (tn / (tn + fp)).fillna(0)
    df["tnr"] = tnr
    precision = (tp / (tp + fp)).fillna(0)
    df["precision"] = precision
    fdr = (fp / (tp + fp)).fillna(0)
    df["fdr"] = fdr
    acc = ((tp + tn) / (tp + tn + fp + fn)).fillna(0)
    df["accuracy"] = acc

    df["youden"] = tpr + tnr - 1
    df["f1_score"] = 2 * (tpr * precision / (tpr + precision)).fillna(0)

    df["mcc"] = (((tp * tn) - (fp * fn)) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))).fillna(0)

    return df


def complete_results(results_df):
    """


    Parameters
    ----------
    results_df

    Returns
    -------

    """

    # Find relations between numerical values and compsition/inrease vectors
    b = []
    for y1_0 in [200, 400, 600, 800, 1000]:
        b_i = np.round(gen.counts_from_first(y1_0, 5000, 5), 3)
        b.append(np.round(np.log(b_i / 5000), 2))

    b_counts = dict(zip([b_i[0] for b_i in b], [200, 400, 600, 800, 1000]))

    b2 = []
    for y1_0 in [200, 400, 600, 800, 1000]:
        b_i = np.round(gen.counts_from_first(y1_0, 5000, 5), 3)
        b2.append(b_i)

    b_w_dict = {}
    w_rel_dict = {}
    i = 0
    for b_i in b2:
        b_t = np.round(np.log(b_i / 5000), 3)
        w_d = {}
        for change in [b_i[0] / 3, b_i[0] / 2, b_i[0], b_i[0] * 2, b_i[0] * 3]:
            _, w = gen.b_w_from_abs_change(b_i, change, 5000)
            w_0 = np.round(w[0], 3)
            w_d[w_0] = change
            rel_change = np.round(change / b_i[0], 2)
            w_rel_dict[w_0] = rel_change
        b_w_dict[b_t[0]] = w_d
        i += 1

    results_df["n_controls"] = [x[0] for x in results_df["n_samples"].tolist()]
    results_df["n_cases"] = [x[1] for x in results_df["n_samples"].tolist()]
    results_df["n_total"] = results_df["n_total"].astype("float")
    results_df["w"] = [x[0][0] for x in results_df["w_true"]]
    results_df["b_0"] = [x[0] for x in results_df["b_true"]]
    results_df["b_count"] = [b_counts[np.round(x, 2)] for x in results_df["b_0"]]

    bs = results_df["b_0"].tolist()
    ws = results_df["w"].tolist()
    increases = []
    rel_changes = []
    for i in range(len(bs)):
        increases.append(b_w_dict[bs[i]][ws[i]])
        rel_changes.append(w_rel_dict[ws[i]])
    results_df["num_increase"] = increases
    results_df["rel_increase"] = rel_changes
    results_df["log_fold_increase"] = np.round(np.log2((results_df["num_increase"] +
                                                        results_df["b_count"]) /
                                                       results_df["b_count"]), 2)

    param_cols = ["b_count", "num_increase", "n_controls", "n_cases", "log_fold_increase"]
    metrics = ["tpr", "tnr", "precision", "accuracy", "youden", "f1_score", "mcc"]
    results_df = results_df.loc[:, param_cols + metrics + ["model"]]

    results_df = results_df.rename(columns={"b_count": "Base", "num_increase": "Increase", "model": "Model",
                                   "log_fold_increase": "log-fold increase"})

    return results_df


def execute_on_server(bash_loction, bash_name, script_location, arguments,
                      python_path="/home/icb/johannes.ostner/anaconda3/bin/python"):
    """
    Script to make a bash script that pushes a job to the ICB CPU servers where another python script is executed.

    Parameters
    ----------
    bash_loction: str
        path to the folder where bash file, out and error files are written
    bash_name: str
        name of your job. The bash file, out and error files will have this name (with respective suffixes)
    script_location: str
        path to the python script
    arguments: list
        list of arguments that is passed to the script

    Returns
    -------

    """

    # Build bash file
    bash_file = bash_loction + bash_name + "_script.sh"
    with open(bash_file, "w") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines(f"#SBATCH -o {bash_loction}{bash_name}_out.o\n")
        fh.writelines(f"#SBATCH -e {bash_loction}{bash_name}_error.e\n")
        fh.writelines("#SBATCH -p cpu_p\n")
        fh.writelines("#SBATCH --constraint='avx'\n")
        fh.writelines("#SBATCH -c 1\n")
        fh.writelines("#SBATCH --mem=32000\n")
        fh.writelines("#SBATCH --nice=10000\n")
        fh.writelines("#SBATCH -t 3-00:00:00\n")

        execute_line = f"{python_path} {script_location} "
        for arg in arguments:
            execute_line = execute_line + str(arg).replace(" ", "") + " "
        fh.writelines(execute_line)

    # Run the bash file you just generated
    os.system(f"sbatch {bash_file}")
