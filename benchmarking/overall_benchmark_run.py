# Run overall benchmark for scCODA paper
import sys
import os

# only need if on server:
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/scCODA/')

# if running on server:
# import benchmarking.benchmark_utils as util
# else:
import benchmarking.benchmark_utils as util


dataset_path = os.path.realpath("../../data/overall_benchmark/generated_datasets_overall_benchmark/")
dataset_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_data/overall_benchmark/generated_datasets_005/')

save_path = os.path.realpath("../../data/benchmark_results/overall_benchmark_results/")
save_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_results/overall_benchmark_005/')

models = ["scCODA"]

# Dont try running this at home! There are 150.000 datasets in this benchmark!
# util.benchmark(dataset_path_server, save_path_server, models, "overall", server=True, keep_sccoda_results=True)
util.benchmark(dataset_path, save_path, models, "overall", server=False, keep_sccoda_results=True)

