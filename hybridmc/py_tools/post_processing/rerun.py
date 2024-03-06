#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# rerun.py reruns the transitions in the rerun dir whose trajectories are
# shorter than the outer mean first passage times
#
# example of how to run (must be in rerun dir):
# python ../../tools/rerun.py ../../release/hybridmc 2 10 -n

#SBATCH --account=def-jmschofi
#SBATCH --time=1-0
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=256M

import sys

sys.path.insert(0, '/home/rgladkik/scratch/experiments/tools')

import argparse
import csv
import glob
import json
import math
import multiprocessing
import os
import subprocess
import re


def main(args):
    csv_in = 'times.csv'

    config_in = []
    config_out = []
    output_name = []
    t = []
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            config_in.append(str(row[0]))
            config_out.append(str(row[1]))
            t.append(float(row[2]))

            layer = str(row[0]).count('1')
            output_name.append('hybridmc_{}_{}_{}'.format(layer, row[0],
                row[1]))

    bits = []
    for file_path in glob.glob('hybridmc_*.h5'):
        bits_from_hdf5 = re.search(r'_([01]+)_([01]+)\.h5$', file_path)
        bits_in, bits_out = bits_from_hdf5.groups()
        bits.append((bits_in, bits_out))

    inputs = []
    for i in range(len(config_in)):
        layer = str(config_in[i]).count('1')

        if layer == 0:
            input_hdf5 = None
        else:
            for (bits_in, bits_out) in bits:
                if bits_out == config_in[i]:
                    input_hdf5 = 'hybridmc_{}_{}_{}.h5'.format(layer-1,
                            bits_in, bits_out)
                    break
            else:
                raise Exception('could not find matching bit pattern for {}, {}'.format(config_in[i], layer))

        inputs.append([i, args.exe, args.seed_increment, args.nsteps,
            input_hdf5, t[i], args.dry_run, output_name[i]])

    with multiprocessing.Pool() as pool:
        pool.starmap(run_layer, inputs)


def run_layer(i, exe, seed_increment, new_nsteps, input_hdf5, t_i, dry_run,
        output_name):
    json_name = '{}.json'.format(output_name)
    hdf5_name = '{}.h5'.format(output_name)

    with open(json_name, 'r') as input_json:
        common_data = json.load(input_json)

    # increase trajectory time such that it is longer than the outer mfpt
    trajectory_t = math.ceil(t_i / 10) * 10 + 10
    common_data['nsteps'] = new_nsteps
    common_data['del_t'] = trajectory_t / new_nsteps

    run_sim(exe, common_data, input_hdf5, hdf5_name, output_name, i,
            seed_increment, dry_run)


def run_sim(exe, data, input_hdf5, hdf5_name, output_name, count,
        seed_increment, dry_run):
    data['seeds'] = [count, seed_increment]

    if os.path.exists(hdf5_name):
        return

    json_name = '{}.json'.format(output_name)
    log_name = '{}.log'.format(output_name)

    # for layer = 1 or greater,
    command = [exe, json_name, hdf5_name]
    if not input_hdf5 is None:
        command += ['--input-file', input_hdf5]

    print(command)
    sys.stdout.flush()

    if dry_run:
        print('would run: {}'.format(command))
        return

    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                stderr=subprocess.STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('seed_increment', type=int,
            help='increment seed with every run')
    parser.add_argument('nsteps', type=int, help='new number of steps')
    parser.add_argument('-n', '--dry-run', action='store_true', default=False,
            help='output commands without running simulations')

    args = parser.parse_args()

    main(args)
