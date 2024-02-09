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
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from py_tools.helpers.run_helpers import init_json
from py_tools.post_processing import avg_s_bias, diff_s_bias, mfpt


def main(args):
    file_name = os.path.basename(args.json)
    dir_name = os.path.splitext(file_name)[0]
    tmp_dir_name = f'{dir_name}.tmp'

    # Check if target directory or tmp already exists. Modify name if needed
    check_dir(args, dir_name)
    check_dir(args, tmp_dir_name)

    # create the temporary directory to run the simulations
    if not os.path.isdir(tmp_dir_name):
        os.mkdir(tmp_dir_name)

    # move into the temporary directory
    os.chdir(tmp_dir_name)

    # Change directory name to suit the new temp working directory.
    # if abs path given not needed


    # add ../ to the path to indicate its use from a directory one more level down.
    (json, exe) = ("../" + path for path in (args.json, args.exe))

    if args.abspath:
        args.json = json
    else:
        args.json, args.exe = json, exe

    # Create dictionary that will have arguments passed to init_json
    init_json_args = {"json": args.json, "seed_increment": 1, "exe": args.exe}

    nproc = os.cpu_count()
    if os.getenv('SLURM_CPUS_PER_TASK'):
        nproc = os.getenv('SLURM_CPUS_PER_TASK')

    init_json_args["nproc"] = nproc

    # run the simulations for the layers
    init_json(init_json_args)

    # Calculate diff sbias, avg sbias and mfpt using simulation results
    post_processing()

    # Move up from the directory with simulation results
    os.chdir("../")

    # Rename the directory -- remove the .tmp tag to show that this simulation has run completely with success
    os.rename(src=tmp_dir_name, dst=dir_name)


def check_dir(args, dir):
    if os.path.isdir(dir):

        print(f'{dir} already exists; saved as old version with given version code or next available code')

        while os.path.exists(f"{dir}_{args.old_version}"):
            args.old_version += 1

        os.rename(src=dir, dst=f"{dir}_{args.old_version}")


def post_processing():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias(out_csv='diff_s_bias.csv')
    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv")
    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt()
    # put together the mfpts in one file
    mfpt.compile_mfpts()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='test_st.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../release/hybridmc")
    parser.add_argument('-ov', '--old_version', help='set version code for old structure simulation run if needed',
                        default=1, type=int)
    parser.add_argument('--abspath', help='set version code for old structure simulation run if needed',
                        default=0, type=int)

    args = parser.parse_args()

    main(args)
