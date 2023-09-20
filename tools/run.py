#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# run.py runs all layers of the layer approach to simulate a protein folding,
# with an optional staircase potential
#
# example of how to run:
# python ../tools/run.py ../../examples/test.json ../../release/hybridmc ../../tools/init_json.py -- --bp 4 16 --rc 3.0 2.2 1.5
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
from helpers.run_helpers import init_json
from post_processing import diff_s_bias, avg_s_bias


def main(args):
    file_name = os.path.basename(args.json)
    dir_name = os.path.splitext(file_name)[0]

    if os.path.isdir(dir_name):
        print(f'{dir_name} already exists; saved as old version with given version code')
        os.rename(src=dir_name, dst=f"{dir_name}_{args.old_version}")

    if os.path.isdir(f"{dir_name}.tmp"):
        os.remove(os.path.realpath(f"{dir_name}.tmp"))

    tmp_dir_name = f'{dir_name}.tmp'

    init_json_args = {"json": args.json, "seed_increment": 1, "exe": args.exe}

    nproc = os.cpu_count()
    if os.getenv('SLURM_CPUS_PER_TASK'):
        nproc = os.getenv('SLURM_CPUS_PER_TASK')

    init_json_args["nproc"] = nproc

    if not os.path.isdir(tmp_dir_name):
        os.mkdir(tmp_dir_name)
    os.chdir(tmp_dir_name)

    init_json(init_json_args)

    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias()

    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv", structure_sim_json=args.json)

    # Move up from the directory with simulation results
    os.chdir("../")

    # Rename the directory -- remove the .tmp tag to show that this simulation has run completely with success
    os.rename(src=tmp_dir_name, dst=dir_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('--old_version', help='set version for old structure simulation run', default='old')

    args = parser.parse_args()

    main(args)
