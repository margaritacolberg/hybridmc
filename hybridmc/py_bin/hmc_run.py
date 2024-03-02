#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# run.py runs all layers of the layer approach to simulate a protein folding,
# with an optional staircase potential
#
# example of how to run:
# python ../tools/run.py --json ../../examples/test.json --exe ../../release/hybridmc
#
# note that run.py creates a new dir which it enters to generate the output
# files, thus, the path of the json input file, init_json.py and the executable
# must be from this newly created dir, and not from the dir from which the
# simulation is initiated

# SBATCH --account=def-jmschofi
# SBATCH --time=3-0
# SBATCH --cpus-per-task=16
# SBATCH --mem-per-cpu=256M

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from hybridmc.py_tools.helpers.run_helpers import init_json
from hybridmc.py_tools.post_processing import avg_s_bias
from hybridmc.py_tools.post_processing import diff_s_bias, constructPaths, mfpt


class RunParams:
    """
    Class to hold the parameters of the hybridmc simulation
    """
    def __init__(self, json, exe, old_version="1", tmp_old_version=""):
        self.json = json
        self.exe = exe
        self.old_version = old_version
        self.tmp_old_version = tmp_old_version


def process(run_params):
    """
    Function to run hybridmc simulation with given RunParams
    """
    file_name = os.path.basename(run_params.json)
    dir_name = os.path.splitext(file_name)[0]
    tmp_dir_name = f'{dir_name}.tmp'

    # Check if target directory or tmp already exists. Modify name if needed
    check_dir(dir_name, old_version_code=run_params.old_version)
    check_dir(tmp_dir_name, old_version_code=run_params.tmp_old_version)

    # Convert the given paths for json and exe to absolute paths if not already absolute
    run_params.json = os.path.realpath(run_params.json)
    run_params.exe = os.path.realpath(run_params.exe)

    # create the temporary directory to run the simulations
    if not os.path.isdir(tmp_dir_name):
        os.mkdir(tmp_dir_name)
    os.chdir(tmp_dir_name)

    # Create dictionary that will have arguments passed to init_json
    init_json_args = {"json": run_params.json, "seed_increment": 1, "exe": run_params.exe}
    nproc = os.cpu_count()
    if os.getenv('SLURM_CPUS_PER_TASK'):
        nproc = int(os.getenv('SLURM_CPUS_PER_TASK'))
    init_json_args["nproc"] = nproc

    # run the simulations for the layers
    init_json(init_json_args)

    # Calculate diff sbias, avg sbias and mfpt using simulation results
    post_processing()

    # Move up from the directory with simulation results
    os.chdir("../../")

    # Rename the directory -- remove the .tmp tag to show that this simulation has run completely with success
    os.rename(src=tmp_dir_name, dst=dir_name)


def check_dir(dir_name, old_version_code=""):
    """
    Checks if the directory exists and do what is specified to do with it:
    - if old_version_code is not given then default to "" then continue last simulation
    - if old_version_code is "delete" restart simulation and delete last results
    - if old_version_code is a number then restart simulation, but tag last results with this number and keep it
    """
    if os.path.isdir(dir_name):

        if old_version_code == "":
            print("Blank old version code given. Using pre-existing tmp directory for simulation")
            return None

        # delete directory if old version code given as delete and return function
        elif old_version_code == "delete":
            print("Removing old tmp directory since old version code supplied as delete")
            return shutil.rmtree(dir_name)

        elif old_version_code.isnumeric():
            old_version_code = int(old_version_code)
            print(f'{dir_name} already exists and int old version code; '
                  f'saved as old version with given version code or next available code')

            # else add available old version code to path name
            while os.path.exists(f"{dir_name}_{args.old_version}"):
                old_version_code += 1

            return os.rename(src=dir_name, dst=f"{dir_name}_{old_version_code}")

        elif os.path.exists(f"{dir_name}{args.old_version}"):
            raise FileExistsError("Tmp file already exists, so does the old version code tmp file")

        else:
            print(f'{dir_name} already exists. Renaming with old version code given')
            return os.rename(src=dir_name, dst=f"{dir_name}{old_version_code}")


def post_processing():
    """
    Post-processing data analysis routine using hybridmc simulation results.
    """
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias(out_csv='diff_s_bias.csv')
    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv", output_csv="avg_s_bias.csv")
    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt()
    # put together the mfpts in one file
    mfpt.compile_mfpts()
    # check if direct and indirect (via paths) diff s bias result consistent
    constructPaths.diff_sbias_state_function_check(diff_s_bias_csv='diff_s_bias.csv', csv_out="diff_check.csv")


def main(run_params):
    """
    Run hybridmc from command line arguments
    """
    process(run_params)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='test.json')
    parser.add_argument('--exe', help='hybridmc executable', default="hybridmc.exe")
    parser.add_argument('-ov', '--old_version', help='set version code for old structure simulation run if needed',
                        default="1", type=str)
    parser.add_argument('-tov', '--tmp_old_version',
                        help='set version code for old structure simulation tmp run if needed',
                        default="", type=str)

    args = parser.parse_args()

    # Create an instance of RunParams from argparse arguments and pass to main
    main(run_params=RunParams(
        json=args.json, exe=args.exe,
        old_version=args.old_version, tmp_old_version=args.tmp_old_version
                   ))
