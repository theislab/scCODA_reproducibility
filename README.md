# scCODA reproducibility

This repository contains the code that was used to produce the results and plots for *scCODA: A Bayesian model for compositional single-cell data analysis* (Büttner et al., 2020). The package containing the model can be found [here](https://github.com/theislab/scCODA).

## Simulated data benchmarks

The directory *benchmarking* contains code for the three benchmarks (threshold determination, model comparison and overall benchmark). The data generation setup for all benchmarks is in `generate_data.py`.
Each benchmark can then be executed via the according `<benchmark_xy>_run.py` script. The analysis and plots are produced in a jupyter notebook named `<benchmark_xy>_analysis.py` for each benchmark.

**Note: Since running the benchmarks is very resource intensive, we provided all generated data and benchmark results at TODO!!!!!!!!!!**

To ensure compatibility with other versions of some python packages, we provided the data and results as .csv and .h5 files instead of .pkl files, which were  produced by the benchmark. Conversion functions are in `benchmarking.to_from_pickle.py`. 

The analysis notebooks can be run with the files `benchmark_results` and `sccoda_effects for each benchmark`. For running the benchmarks, please convert the files `generated_data` and `generation_parameters` via the function `benchmarking.to_from_pickle.benchmark_datasets_to_pickle`.

## Applications

All applications of scCODA that were shown in the article can be found under *applications*. For each analysis, a jupyter notebook is provided.
