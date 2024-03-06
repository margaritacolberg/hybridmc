#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# get_rerun_files.py finds the transitions to rerun whose trajectories are
# shorter than the outer mean first passage times; the bit patterns and outer
# mfpt are written to an output file, and the corresponding json and hdf5 files
# are moved to the rerun dir
#
# example of how to run:
# python ../tools/get_rerun_files.py ../examples/test.json rerun

import sys

sys.path.insert(0, '/home/rgladkik/scratch/experiments/tools')
import matrix_element

import argparse
import csv
import json
import glob
import os

from shutil import copyfile


def main(args):
    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    nsteps = data['nsteps']
    del_t = data['del_t']

    os.makedirs(args.dirpath, exist_ok=True)

    bits = []
    config_in = []
    output = []
    for file_path in glob.glob('hybridmc_*.csv'):
        bits_i, bits_j, t_ind, _int_i, _int_j = \
                matrix_element.get_state_data(file_path)
        bits.append((bits_i, bits_j))

        inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)

        # length of trajectory is nsteps * del_t
        if outer_fpt > (nsteps * del_t):
            layer = str(bits_i).count('1')

            json_name = 'hybridmc_{}_{}_{}.json'.format(layer, bits_i, bits_j)
            copyfile(json_name, os.path.join(args.dirpath, json_name))

            config_in.append(bits_i)
            output.append([bits_i, bits_j, outer_fpt])

    get_hdf5_input(config_in, bits, args.dirpath)

    csv_out = os.path.join(args.dirpath, 'times.csv')
    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def get_hdf5_input(config_in, bits, dirpath):
    for i in range(len(config_in)):
        layer = str(config_in[i]).count('1')

        if layer == 0:
            input_hdf5 = None
        else:
            for (bits_i, bits_j) in bits:
                if bits_j == config_in[i]:
                    input_hdf5 = 'hybridmc_{}_{}_{}.h5'.format(layer-1, bits_i,
                            bits_j)

                    for hdf5_file in glob.glob(input_hdf5):
                        copyfile(hdf5_file, os.path.join(dirpath, hdf5_file))

                    break
            else:
                raise Exception('could not find matching bit pattern for {}'.format(config_in[i]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('dirpath', help='directory path to store rerun files')

    args = parser.parse_args()

    main(args)
