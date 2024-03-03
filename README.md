# HybridMC: Hybrid Monte Carlo Protein Folding

[![Anaconda-Server Badge](https://anaconda.org/vignesh229/hybridmc/badges/version.svg)](https://anaconda.org/vignesh229/hybridmc)
[![Anaconda-Server Badge](https://anaconda.org/vignesh229/hybridmc/badges/platforms.svg)](https://anaconda.org/vignesh229/hybridmc)
[![Publish Conda Package on Linux](https://github.com/Vigne-hub/hybridmc/actions/workflows/conda_publish_linux.yml/badge.svg)](https://github.com/Vigne-hub/hybridmc/actions/workflows/conda_publish_linux.yml)

HybridMC simulates the event-driven dynamics of coarse-grained protein folding.
The entropy and mean first passage times (MFPT) of each bond forming or
breaking event are calculated, and the Markov transition rate matrix is
constructed. The average time needed to fold to the protein's native state,
starting from the unfolded state, is evaluated under two conditions.

For more information, see the [publication](https://doi.org/10.1063/5.0098612)
or the [preprint](https://arxiv.org/abs/2205.05799).

## Directories

See top level comment in each file for more details as to what each file does.
Some top level comments also contain examples of how to run each file, if
complicated.

  * `src`: C++ source code for HybridMC program

  * `examples`: JSON inputs which specify protein configuration and Bash
    scripts for running HybridMC

  * `py_tools`: Python module with some utilituy functions to automate running HybridMC and extracting data such as
    the MFPT

  * `py_bin`: Python scripts for running HybridMC

  * `tree`: C++ and Bash scripts for generating tree plots

  * `crambin_s_bias_mfpt`: biased entropy and MPFT results for crambin

## Installing Package

This package can be installed using **`conda install hybridmc -c conda-forge -c vignesh229`** to directly obtain a command called **hmc_run** which can be used simply as **`hmc_run --json INPUT_JSON_FILE`**. This also installs all the python utilities that can be imported from a python shell.

## C++ Program Details for Manual Install
**NOTE:** This is not the recommended approach, proceeed to use this with caution.

### Prerequisite Software

  * C++17 compiler with support for C++20
    [bit library](https://en.cppreference.com/w/cpp/header/bit)

    see bit operations under C++20 [compiler
    support](https://en.cppreference.com/w/cpp/compiler_support/20)

    tested with GCC 9.3.0 on CentOS 7.9

  * Minimum CMake version is 3.13.1

  * [Ninja](https://ninja-build.org/)

  * HDF5 version 1.10 or later

  * Boost

All of these can be obtained through the conda environment yaml file: hybridmc.yaml

You can get these within a conda environment by running:

```
conda install cmake hdf5 ninja cxx-compiler boost -c conda-forge
```

### Commands for Compiling

To compile HybridMC for development and testing:

```
make debug
```

To perform unit test with verbose output:

```
cd debug
ctest -V
```

To compile HybridMC for running simulations:

```
make release
```
To then install the wang landau plugin for python, go to the main hybridmc directory and run

```
pip install .
```

## Python Program Details

### Run Program

Simply run the python executable run.py --json JSON_NAME.json

The json files for running the examples (for example test.json) in the examples directory
have the following parameters:


    "m": 1.0, -- mass of each bead

    "sigma_bb": 1.0, -- 

    "near_min": 1.0, -- shortest bond dist between nearest neighbors

    "near_max": 1.17, -- longest bond dist between nearest neighbors

    "nnear_min": 1.4, -- shortest bond dist between next nearest neighbors

    "nnear_max": 1.67, -- longest bond dist between next nearest neighbors

    "rh": 1.25, -- The minimum distance between nonlocal beads 

    "rc": 1.5, -- the minimum distance between transient beads. Default value top add to list of bonds if not specified
    
    "nonlocal_bonds": [[2, 6], [6, 10, 1.7], [10, 14]], -- list of  nonlocal beads. Could specify rc for each bond
    
    "transient_bonds": [[2, 6], [6, 10], [10, 14]], -- list of transient beads
    
    "permanent_bonds": [], -- list of permanent beads
    
    "config_in": 0, 

    "config_out": 7,

    "nbeads": 20, -- number of beads in the protein

    "tries": 10000, -- number of tries to find a valid configuration

    "length": 18.0, -- length of the periodic box in each dimension

    "ncell": 5, -- number of cells in each dimension

    "nsteps": 8, -- number of steps in each iteration

    "nsteps_eq": 8, -- number of steps in each iteration for equilibration

    "del_t": 15.0, -- length of dynamical trajectory used for sampling

    "nsteps_wl": 1, -- number of steps for wang_landau

    "del_t_wl": 4.0, -- length of dynamical trajectory used for wang_landau

    "gamma": 0.1, -- starting adjustment factor for wang landau

    "gamma_f": 0.001, -- final adjustment factor for wang landau

    "gamma_f_screening": 0.0005, -- wang landau screening adjustment factor

    "seeds": [3, 1], -- seeds for random number generator

    "temp": 1.0, -- temperature

    "mc_moves": 20, -- number of monte carlo moves

    "total_iter": 400, -- number of state counts per iteration used to calculate the sbias. Increasing this imporves sbias calculation accuracy
    
    "total_iter_eq": 100, -- total number of iterations for equilibration

    "pos_scale": 0.5,

    "neg_scale": 0.1,

    "sig_level": 0.05, -- significance level for the g-test

    "max_nbonds": 3, -- maximum number of bonds for the protein

    "max_g_test_count": 10000, -- maximum number of g-tests to run

    "flip_req": 1.0, -- fraction of bonds that must be flipped to consider a configuration valid

    "fail_max": 5, -- maximum number of times to fail to find a valid configuration before increasing the numebr of iterations

    "WL_sbias": 6.0, -- the sbias value for a wang landau simulation result beyond which a configuration is deemed to need staircasing 

    "req_dists": 50000, -- minimum distances used for MFPT. Enforced in the hybridmc c++ source: it is a convergence criterion.

    "rc_target_min_percentile": 0.025, -- Determines the lcoation of rc for the stair. Look at 2.5 percentile rc value from the rc explored
                                            during simulation.

    "D_mean": 0.046, -- Estimated diffusion coefficent mean value. Used for rate calculation.
    
    "D_std_perc": 1, -- Assumed % error in D_mean.

    "eps": 3, -- The energy of the bond in kT units. Used to calculate reates and probability.

    "useEnsemble": false, -- Experimental feature switch: swapping states during initialization
    
    "ensembleSize": 100, -- If use enemble how many states to choose from
    
    "ensemble_write_step": 10 -- If use ensemble, steps after which to swap in new state

}

### Calculate Entropy

To get biased entropy of each state manually, use `diff_s_bias.py`, followed by
`avg_s_bias.py`. To get biased entropy of each state which includes the
staircase potential, use `diff_s_bias_stair.py`, followed by `avg_s_bias.py`.

But this is done automatically when run.py is executed.

To find which transition gives smallest entropy difference, use
`min_diff_s_bias.py`. To get percent error for the entropy of one transition,
use `get_error_s_bias.py`.

In addition, an output for a state function check if provided for the diff_s_bias (diff_check.csv)
with information about if the direct path value is significantly different
than an average value using multiple paths for the same transition.

### Calculate MFPT

To calculate MFPT for a single transition manually, use `mfpt.py`. To calculate MFPT for
all transitions in current dir, use `run_mfpt.py`. To combine the MFPT for each
step of staircase into a single MFPT for the transition, use
`mfpt_for_stair.py`.

Again, as with the entropy, this is done automatically when run.py is executed.

### Calculate Energy and Folding Time

To calculate energy of each intermediate state and average time to fold to the
native state:

 1. use `fixed_e.py` if all nonlocal bonds have the same energy, and
    temperature varies with each repeat of the simulation, or

 2. use `min_e.py` to hold probability of native state high relative to the
    unfolded state.

### Tree Plots

To output data in the correct format for the tree plots, run `fixed_e_tree.py`
or `min_e_tree.py` to create a .dat and .mtx file. Then, to generate the tree
plots for crambin, run `crambin_fixed_e.sh` or `crambin_min_e.sh`, if the
executables disconnectionDPS and constructPmatrix exist.

To create the disconnectionDPS executable, the disconnectionDPS package is
needed, which can be downloaded
[here](https://www-wales.ch.cam.ac.uk/software.html). The download contains
many subdirectories, but the one needed for the tree plots is called
DISCONNECT. Inside DISCONNECT/source/, to compile disconnectionDPS, run

```
make F90=gfortran
```

To compile constructPmatrix, run

```
g++ -O2 -o constructPmatrix constructPmatrix.cc $(gsl-config --cflags --libs)
```

### Debug

To debug one transition which compiled successfully but crashed once it started
running, use `testing.py`.
