#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# get_error_s_bias.py runs the dynamics of the same transition repeatedly to
# obtain a set of biased entropies for which the percent error is evaluated
#
# example of how to run:
# python ../tools/get_error_s_bias.py hybridmc_1_001_011.json ../release/hybridmc --hdf5 hybridmc_0_000_001.h5 10

#SBATCH --account=def-jmschofi
#SBATCH --time=1-0
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=256M

import argparse
import csv
import h5py
import json
import multiprocessing
import numpy as np
import os
import subprocess
import sys


class OutputEntropy(object):
    def __init__(self, data, output_name, exe, hdf5, seed_increment):
        self.data = data
        self.output_name = output_name
        self.exe = exe
        self.hdf5 = hdf5
        self.seed_increment = seed_increment

    def __call__(self, i):
        return get_s_bias(self.data, self.output_name, self.exe, self.hdf5,
                self.seed_increment, i)


def main(args):
    input_name = os.path.basename(args.json)
    output_name = os.path.splitext(input_name)[0]

    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    output_entropy = OutputEntropy(data, output_name, args.exe, args.hdf5,
            args.seed_increment)
    with multiprocessing.Pool() as pool:
        s_bias = pool.map(output_entropy, range(args.nboot))

    s_bias_mean = np.mean(s_bias)
    print('mean of entropy samples:', s_bias_mean)

    s_bias_var = np.var(s_bias)
    print('var of entropy samples:', s_bias_var)

    rel_e = s_bias_var / s_bias_mean
    percent_rel_e = rel_e * 100
    print('percent error:', percent_rel_e)

    csv_name = '{}_s_bias_error.csv'.format(output_name)
    output = [[s_bias_mean, s_bias_var, percent_rel_e]]
    with open(csv_name, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def get_s_bias(data, output_name, exe, hdf5_in, seed_increment, i):
    hdf5_name = '{}_{}.h5'.format(output_name, i)
    json_name = '{}_{}.json'.format(output_name, i)

    if not os.path.exists(hdf5_name):
        # each json file must have a unique seed
        data['seeds'] = [i, seed_increment]

        with open(json_name, 'w') as input_json:
            json.dump(data, input_json)

        command = [exe, json_name, hdf5_name]
        if not hdf5_in is None:
            command += ['--input-file', hdf5_in]

        print(command)
        sys.stdout.flush()

        log_name = '{}_{}.log'.format(output_name, i)
        with open(log_name, 'w') as output_log:
            subprocess.run(command, check=True, stdout=output_log,
                    stderr=subprocess.STDOUT)

    with h5py.File(hdf5_name, 'r') as f:
        s_bias = (f['s_bias'][0])

    return s_bias


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('--hdf5', help='hdf5 input file')
    parser.add_argument('nboot', type=int, help='number of bootstraps')
    parser.add_argument('--seed_increment', type=int, default=2,
            help='increment seed with every run')

    args = parser.parse_args()

    main(args)
