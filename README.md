# scCODA reproducibility

This repository contains the code that was used to produce the results and plots for *scCODA: A Bayesian model for compositional single-cell data analysis* (Büttner et al., 2020). The package containing the model can be found [here](https://github.com/theislab/scCODA).

## Simulated data benchmarks

The directory *benchmarking* contains code for the three benchmarks (threshold determination, model comparison and overall benchmark). The data generation setup for all benchmarks is in `generate_data.py`.
Each benchmark can then be executed via the according `<benchmark_xy>_run.py` script. The analysis and plots are produced in a jupyter notebook named `<benchmark_xy>_analysis.py` for each benchmark.

**Note: Since running the benchmarks is very resource intensive, we provided all generated data and benchmark results at TODO!!!!!!!!!!**

To ensure compatibility with other versions of some python packages, we provided the data and results as .csv and .h5 files instead of .pkl files, which are produced by the code. Conversion functions are in `benchmarking.to_from_pickle.py`. 

The analysis notebooks can be run with the files `benchmark_results` and `sccoda_effects for each benchmark`. For running the benchmarks, please convert the files `generated_data` and `generation_parameters` via the function `benchmarking.to_from_pickle.benchmark_datasets_to_pickle`:

### Example workflow for reproducing a benchmark

1. Navigate to the parent directory of `scCODA_reproducibility`

2.  a) Download the benchmark data and results from **TODO!!!!** and unpack them here

    **or**
    
    b) If you want to re-run a benchmark, create the following directory structure:
    
    ```
    <parent_directory>
    |_  <scCODA_reproducibility>
    |_  data
        |_  overall_benchmark
            |_  data_overall_benchmark
            |_  overall_benchmark_results
            |_  generated_datasets_overall_benchmark
            |_  overall_benchmark_plots
        |_  model_comparison_benchmark
            |_  data_model_comparison
            |_  model_comparison_results
            |_  generated_datasets_model_comparison
            |_  model_comparison_plots
        |_  threshold_determination_benchmark
            |_  data_threshold_determination
            |_  threshold_determination_results
            |_  generated_datasets_threshold_determination
            |_  threshold_determination_plots
    ```
    
    Then run `benchmarking.generate_data` to generate all benchmark data. 
    The `generated_datasets_<benchmark_name>` directories should now be filled with files.
    
    Then, execute the according script `benchmarks.<benchmark_name>_run`.
    This produces results in the directory `<benchmark_name>_results`.

    Once you finished re-running all the benchmarks you want, run `benchmarking.to_from_pickle` to convert the data to standardized file types.
    
3. To re-do the benchmark analysis, run the according jupyter notebook `benchmarking.<benchmark_name>_analysis`.


## Applications

All applications of scCODA that were shown in the article can be found under *applications*. For each analysis, a jupyter notebook is provided.
