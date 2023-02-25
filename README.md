# HybridMC: Hybrid Monte Carlo Protein Folding

[![build workflow status](https://github.com/margaritacolberg/hybridmc/actions/workflows/build.yml/badge.svg)](https://github.com/margaritacolberg/hybridmc/actions/workflows/build.yml?query=branch:main)
[![format workflow status](https://github.com/margaritacolberg/hybridmc/actions/workflows/format.yml/badge.svg)](https://github.com/margaritacolberg/hybridmc/actions/workflows/format.yml?query=branch:main)

HybridMC simulates the event-driven dynamics of coarse-grained protein folding.
The entropy and mean first passage times (MFPT) of each bond forming or
breaking event are calculated, and the Markov transition rate matrix is
constructed. The average time needed to fold to the protein's native state,
starting from the unfolded state, is evaluated under two conditions.

For more information, see [preprint](https://arxiv.org/abs/2205.05799).

## Directories

See top level comment in each file for more details as to what each file does.
Some top level comments also contain examples of how to run each file, if
complicated.

  * `src`: C++ source code for HybridMC program

  * `examples`: JSON inputs which specify protein configuration and Bash
    scripts for running HybridMC

  * `tools`: Python scripts for running HybridMC and extracting data such as
    the MFPT

  * `tree`: C++ and Bash scripts for generating tree plots

  * `crambin_s_bias_fpt`: biased entropy and MPFT results for crambin

## C++ Program Details

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

## Python Program Details

### Run Program

For any protein with or without staircase potentials, use `run.py`; for
non-staircase crambin folding, use `run.py` or `crambin.sh`.

### Calculate Entropy

To get biased entropy of each state, use `diff_s_bias.py`, followed by
`avg_s_bias.py`. To get biased entropy of each state which includes the
staircase potential, use `diff_s_bias_stair.py`, followed by `avg_s_bias.py`.

To find which transition gives smallest entropy difference, use
`min_diff_s_bias.py`. To get percent error for the entropy of one transition,
use `get_error_s_bias.py`.

### Calculate MFPT

To calculate MFPT for a single transition, use `mfpt.py`. To calculate MFPT for
all transitions in current dir, use `run_mfpt.py`. To combine the MFPT for each
step of staircase into a single MFPT for the transition, use
`mfpt_for_stair.py`.

### Rerun Files

To move JSON and HDF5 files needed for rerunning transitions to rerun dir, use
`get_rerun_files.py`. To rerun files that were moved to rerun dir, use
`rerun.py`.

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
