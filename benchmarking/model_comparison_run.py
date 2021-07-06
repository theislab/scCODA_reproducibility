# Run model comparison benchmark for scCODA paper
import sys
import os

# only need if on server:
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/')

import benchmarking.benchmark_utils as util

dataset_path = os.path.realpath("../sccoda_benchmark_data/model_comparison_2/generated_datasets_model_comparison_2/") + "/"
dataset_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_data/generated_datasets_model_comparison_2/') + "/"

save_path = os.path.realpath("../sccoda_benchmark_data/model_comparison_2/model_comparison_2_results/") + "/"
save_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_results/model_comparison_2/') + "/"

# Use all 10 models for comparison
models_server = ["scCODA", "scdc", "simple_dm"]
models_local = ["alr_ttest", "ALDEx2_alr", "alr_wilcoxon", "Haber", "ttest", "ancom", "dirichreg", "ANCOMBC_holm", "BetaBinomial"]

util.benchmark(dataset_path, save_path, models_local, "comp_2", server=False, keep_results=True)
# util.benchmark(dataset_path_server, save_path_server, models_server, "comp_2", server=True, keep_results=True)
