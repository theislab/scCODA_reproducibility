# Run threshold benchmark for scCODA paper
import sys
import os

# only need if on server:
# sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/scCODA/')

# if running on server:
# import benchmark_utils as util
# else:
import benchmarking.benchmark_utils as util


dataset_path = os.path.realpath("../../data/threshold_determination/generated_datasets_threshold_determination/")
dataset_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_data/threshold_determination/generated_datasets_005_balanced/')

save_path = os.path.realpath("../../data/threshold_determination/threshold_determination_results/")
save_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_results/threshold_determination_005_balanced/')

models = ["scCODA"]

util.benchmark(dataset_path, save_path, models, "threshold", server=False, keep_sccoda_results=True)
# util.benchmark(dataset_path_server, save_path_server, models, "threshold", server=True, keep_sccoda_results=True)
