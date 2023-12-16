## Early warning signs of critical transitions - The alpha-stable case (Layritz et al, 2023)

This code reproduces the data and analysis for https://arxiv.org/pdf/2311.16350.pdf

The simulation code was run on the CoolMUC2 Cluster of the Leibniz Supercomputing Center. The output of these simulations are found in `data`.

If you just want to test the simulations locally, `test_run.py` will generate one iteration of both systems for both the linear and non-linear case. The non-linear case is commented out by default as this takes some time to run.

For this run `python3 code/run_test.py`. The output will appear in `data/new`

Running `figures.R` will produce all figures. The figures will be found in `reports/paper`.

Conda environment

```
conda create --name alphastableEWS python=3.9
source activate alphastableEWS
conda install -c conda-forge numpy=1.23.1 scipy=1.9.0 pandas=1.4.3 
```


