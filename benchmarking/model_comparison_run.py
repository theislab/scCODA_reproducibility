# Run model comparison benchmark for scCODA paper
import sys
import os

# only need if on server:
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/scCODA/')

# if running on server:
# import benchmarking.benchmark_utils as util
# else:
import benchmarking.benchmark_utils as util

dataset_path = os.path.realpath("../../data/model_comparison/generated_datasets_model_comparison/")
dataset_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_data/model_comparison/scdc_missing_datasets/')

save_path = os.path.realpath("../../data/model_comparison/model_comparison_results/")
save_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_results/scdc_completion/')

# Use all 10 models for comparison
models = ["simple_dm", "alr_ttest", "ALDEx2_alr", "alr_wilcoxon", "Haber", "ttest", "scCODA", "ancom", "scdc", "dirichreg"]

util.benchmark(dataset_path, save_path, models, "comp", server=False)
# util.benchmark(dataset_path_server, save_path_server, models, "comp", server=True)
