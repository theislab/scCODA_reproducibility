# Run model sparse effect test
import sys
import os

# only need if on server:
sys.path.insert(0, '/home/icb/johannes.ostner/compositional_diff/')

import benchmarking.benchmark_utils as util

dataset_path = os.path.realpath("../../../other_data/sparsity_test/generated_datasets_model_sparsity_test/") + "/"
dataset_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_data/generated_datasets_model_sparsity_test/') + "/"

save_path = os.path.realpath("../../../other_data/sparsity_test/sparsity_test_results/") + "/"
save_path_server = os.path.realpath('/home/icb/johannes.ostner/compositional_diff/benchmark_results/model_sparsity_2/') + "/"

models_server = ["scCODA"]

# util.benchmark(dataset_path, save_path, models_local, "comp_2", server=False, keep_results=True)
util.benchmark(dataset_path_server, save_path_server, models_server, "sparsity", server=True, keep_results=True)
