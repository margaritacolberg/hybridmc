#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# avg_s_bias.py adds the entropy differences of all transitions from the native
# state to the current state, and this sum becomes the entropy of the current
# state (the entropy of the native state is 0); whenever two (or more)
# transitions lead to the same current state, an average of the entropy values
# is taken
#
# example of how to run:
# python ../tools/avg_s_bias.py ../examples/crambin.json diff_s_bias.csv

import argparse
import csv
import json
import numpy as np


def main(args):
    with open(args.csv, 'r') as input_csv:
        reader = csv.reader(input_csv, delimiter=',')
        data_csv = list(reader)

    diff_s_bias = dict()
    for i in range(len(data_csv)):
        diff_s_bias[(data_csv[i][0], data_csv[i][1])] = float(data_csv[i][2])

    with open(args.json, 'r') as input_json:
        data_json = json.load(input_json)

    nonlocal_bonds = data_json['nonlocal_bonds']
    nbonds = len(nonlocal_bonds)

    # initially fully bonded
    work_list = [[True]*nbonds]

    bonded_config = format_bits(work_list[0])

    s_bias = dict()
    s_bias[bonded_config] = 0

    avg_s_bias = dict()

    while work_list:
        bits_in = work_list.pop(0)
        config_in = format_bits(bits_in)
        mean_s_bias_in = np.mean(s_bias[config_in])
        avg_s_bias[config_in] = mean_s_bias_in

        for i in range(nbonds):
            # if two beads are not bonded,
            if not bits_in[i]:
                # skip; do not form the bond
                continue

            # copy of initial state for every bonded pair of beads
            bits_out = bits_in.copy()
            # flip bit to break bond
            bits_out[i] = False

            config_out = format_bits(bits_out)
            diff_s_bias_out = diff_s_bias[(config_out, config_in)]

            s_bias_out = s_bias.setdefault(config_out, [])
            s_bias_out.append(diff_s_bias_out + mean_s_bias_in)

            # check if bond pattern already exists inside worklist;
            # only append unique bond patterns to worklist
            if not bits_out in work_list:
                work_list.append(bits_out)

    with open('avg_s_bias.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(sorted(avg_s_bias.items()))


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('csv', help='csv input file')

    args = parser.parse_args()

    main(args)
