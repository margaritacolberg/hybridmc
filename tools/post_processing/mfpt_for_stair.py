#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# mfpt_for_stair.py determines the mean first passage times for the transitions
# with staircase potentials; the inner fpt (the first time at which two beads
# at a distance r > rc diffuse to rc) is taken to be the inner fpt of the step
# whose rc is smallest, and the outer fpt (with r > rc) is taken from the
# previous layer whose final state is equal to the initial state of the current
# layer
#
# example of how to run:
# python ../tools/mfpt_for_stair.py hybridmc_0_000_001.json ../release/hybridmc stair_csv
#
# note: mfpt.py must be run for all transitions first, before mfpt_for_stair.py

#SBATCH --account=def-jmschofi
#SBATCH --time=3-0
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=256M

import sys

sys.path.insert(0, '/home/rgladkik/scratch/experiments/tools')
import matrix_element
import mfpt
import argparse
import csv
import glob
import json
import numpy as np
import os
import re
import shutil
import subprocess


def main(args):
    os.makedirs(args.dirname, exist_ok=True)

    files_to_move = []
    stair_in = {}
    stair_out = {}
    for file_path in glob.glob('hybridmc_*_*_*_*.csv'):
        print(file_path)
        bits_from_csv = \
                re.search(r'_(?P<bits_i>[01]+)_(?P<bits_j>[01]+)(?:_(?P<stair>[0-9]+\.[0-9]+))?\.csv$',
                file_path)
        bits_i, bits_j, rc = bits_from_csv.groups()

        if bits_from_csv is None:
            continue

        stair_in.setdefault((bits_i, bits_j), None)
        stair_out.setdefault((bits_i, bits_j), [])

        # excludes file for last step in staircase
        files_to_move.append(file_path)

    # last step in staircase
    for keys, t in stair_in.items():
        last_stair = 'hybridmc_{}_{}_{}.csv'.format(str(keys[0]).count('1'),
                keys[0], keys[1])

        t_ind = [j for j in range(len(keys[0])) if keys[0][j] != keys[1][j]][0]
        inner_fpt, _outer_fpt = matrix_element.get_fpt(last_stair, t_ind)

        # set the inner fpt for the staircase to be the inner fpt of the step
        # whose rc is the smallest
        stair_in[keys] = inner_fpt

        files_to_move.append(last_stair)

    # after inner fpt are extracted, the csv files for each step are no longer
    # needed so they are moved to dirname
    for i in range(len(files_to_move)):
        shutil.move(files_to_move[i], args.dirname)

    bits = []
    config_in = []
    config_out = []
    for file_path in glob.glob('hybridmc_*.csv'):
        bits_i, bits_j, _t_ind, _int_i, _int_j = \
                matrix_element.get_state_data(file_path)
        bits.append((bits_i, bits_j))

    # set the outer fpt for the staircase as the outer fpt from the previous
    # layer whose final state is the initial state of the current layer
    for keys, t in stair_out.items():
        for (bits_i, bits_j) in bits:
            if bits_j == keys[0]:
                file_path = \
                'hybridmc_{}_{}_{}.csv'.format(str(bits_i).count('1'), bits_i,
                        bits_j)

                t_ind = [j for j in range(len(keys[0])) if keys[0][j] !=
                        keys[1][j]][0]
                inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)

                stair_out[keys].append(outer_fpt)

    output = []
    for keys, t in stair_out.items():
        if t:
            stair_out[keys] = np.max(stair_out[keys])
            output.append([keys[0], keys[1], stair_in[keys], stair_out[keys]])
        else:
            # for 0 layer, repeat the transition but without the staircase
            # potential (in the case where a bond occurs between beads located
            # on opposite ends of the chain, convergence will be slow, but in
            # the case where permanent bonds in the initial state may introduce
            # constraints along the chain, convergence will be fast since layer
            # 0 has no permanent bonds)
            with open(args.json, 'r') as input_json:
                data = json.load(input_json)

            del data['stair_bonds']
            del data['stair']

            name = os.path.splitext(args.json)[0]
            json_name = name + '_new.json'
            hdf5_name = name + '_new.h5'
            log_name = name + '_new.log'
            csv_name = name + '.csv'

            with open(json_name, 'w') as output_json:
                json.dump(data, output_json)

            command = [args.exe, json_name, hdf5_name]

            print(command)
            sys.stdout.flush()

            with open(log_name, 'w') as output_log:
                subprocess.run(command, check=True, stdout=output_log,
                        stderr=subprocess.STDOUT)

            mfpt.fpt_write(json_name, hdf5_name, 0, csv_name)

    for i in range(len(output)):
        layer = str(output[i][0]).count('1')
        csv_out = 'hybridmc_{}_{}_{}.csv'.format(layer, output[i][0],
                output[i][1])
        t_ind = [j for j in range(len(output[i][0])) if output[i][0][j] !=
                output[i][1][j]][0]

        with open(csv_out, 'w') as output_csv:
            writer = csv.writer(output_csv)
            writer.writerow([t_ind, output[i][2], output[i][3]])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('dirname',
            help='dir name to store staircase csv files')

    args = parser.parse_args()

    main(args)
