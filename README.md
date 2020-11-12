# scCODA reproducibility

This repository contains the code that was used to produce the results and plots for *scCODA: A Bayesian model for compositional single-cell data analysis* (Ostner et al., 2020). The package containing the model can be found [here](https://github.com/theislab/scCODA).

The directory *benchmarking* contains code for the three benchmarks (threshold determination, model comparison and overall benchmark). The data generation setup for all benchmarks is in `generate_data.py`. Each benchmark can then be run with the according `<benchmark_xy>_run.py` script. The analysis and plots are produced in a jupyter notebook named `<bechmark_xy>_analysis.py` for each benchmark.

**Note: Since running the benchmarks is very resource intensive, we provided all generated data and benchmark results at TODO!!!!!!!!!!**

All applications of scCODA that were shown in the article can be found under *applications*. For each application, a jupyter notebook is provided.
