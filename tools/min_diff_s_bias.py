#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# min_diff_s_bias.py finds the 5 smallest entropy differences for each bond
# that is turned on

import argparse
import csv
import json
import numpy as np


def main(args):
    bits_in = []
    bits_out = []
    s_bias = []

    with open(args.json, 'r') as input_json:
        json_data = json.load(input_json)

    nonlocal_bonds = json_data['nonlocal_bonds']
    nbonds = len(nonlocal_bonds)

    with open(args.csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            bits_in.append(row[0])
            bits_out.append(row[1])
            s_bias.append(float(row[2]))

    output = []
    for i in range(nbonds):
        bits_in_i = []
        bits_out_i = []
        s_bias_i = []
        for j in range(len(s_bias)):
            t_ind = [k for k in range(len(bits_in[j])) if bits_in[j][k] !=
                    bits_out[j][k]][0]

            if i == t_ind:
                bits_in_i.append(bits_in[j])
                bits_out_i.append(bits_out[j])
                s_bias_i.append(s_bias[j])

        # find transition with smallest entropy difference
        min_output = sorted(zip(bits_in_i, bits_out_i, s_bias_i), key=lambda t:
                t[2])[:5]

        for j in range(len(min_output)):
            output.append(np.array(min_output[j]))

    with open(args.csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('csv_in', help='diff_s_bias.csv file')
    parser.add_argument('csv_out', help='csv output file')

    args = parser.parse_args()

    main(args)
