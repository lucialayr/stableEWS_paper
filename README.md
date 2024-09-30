## A Generalized Approach to Early Warning Signs of Critical Transitions in Dynamical Systems driven by α-Stable Noise

This code reproduces the data and analysis for the paper 'A Generalized Approach to Early Warning Signs of Critical Transitions in Dynamical Systems driven by α-Stable Noise' currently under considerations at Theoretical Ecology

The simulation code was optimized for efficient execution on parallelized High-Performance Computing (HPC) infrastructure. The output of these simulations is found in `data`.

If you just want to test the simulations locally, `test_run.py` will generate one iteration of both systems for the linear and non-linear cases. The non-linear case is commented out by default as this takes some time to run.

For this run `python3 code/run_test.py`. The output will appear in `data/new`

Running `figures.R` will produce all figures. The figures will be found in `reports/paper`.

Conda environment

```
conda create --name alphastableEWS python=3.9
source activate alphastableEWS
conda install -c conda-forge numpy=1.23.1 scipy=1.9.0 pandas=1.4.3 
```


