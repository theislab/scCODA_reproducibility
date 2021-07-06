# Only relevant for server execution
import pickle as pkl
import sys
import os

import benchmark_utils as add

dataset_path = sys.argv[1]
save_path = sys.argv[2]
model_name = sys.argv[3]
alpha = float(sys.argv[4])
count = int(sys.argv[5])
if sys.argv[6] == "True":
    keep_results = True
else:
    keep_results = False
print("model name:", model_name)

file_name = os.listdir(dataset_path)[count]

r_home = "/home/icb/johannes.ostner/anaconda3/lib/R"
r_path = r"/home/icb/johannes.ostner/anaconda3/lib/R/bin"

if model_name == "ALDEx2_alr":
    kwargs = {"server": True,
              "method": "we.eBH",
              "mc_samples": 128,
              "denom": [5],
              "alpha": alpha,
              "fdr_correct": False,
              "keep_results": keep_results}

elif model_name == "ALDEx2":
    kwargs = {"server": True,
              "method": "we.eBH",
              "mc_samples": 128,
              "alpha": alpha,
              "fdr_correct": False,
              "keep_results": keep_results}

elif model_name in ["simple_dm", "scCODA"]:
    kwargs = {"num_results": 20000,
              "num_burnin": 5000,
              "alpha": alpha,
              "num_adapt_steps": 4000,
              "keep_results": keep_results}

elif model_name in ["alr_ttest", "alr_wilcoxon"]:
    kwargs = {"reference_index": 4,
              "alpha": alpha,
              "fdr_correct": True,
              "keep_results": keep_results}
elif model_name in ["Haber", "ttest", "clr_ttest", "dirichreg"]:
    kwargs = {"alpha": alpha,
              "fdr_correct": True,
              "keep_results": keep_results}
elif model_name == "scdc":
    kwargs = {
        "r_home": r_home,
        "r_path": r_path,
        "keep_results": keep_results,
        "alpha": alpha
    }
else:
    kwargs = {}

if not dataset_path[-1] == "/":
    dataset_path = dataset_path + "/"

if not save_path[-1] == "/":
    save_path = save_path + "/"

if keep_results:
    results, effects = add.model_on_one_datafile(dataset_path+file_name, model_name, **kwargs)
    results = add.get_scores(results)
    save = {"results": results, "effects": effects}
else:
    results = add.model_on_one_datafile(dataset_path+file_name, model_name, **kwargs)
    results = add.get_scores(results)
    save = results

with open(save_path + model_name + "_results_" + str(count) + ".pkl", "wb") as f:
    pkl.dump(save, f)
