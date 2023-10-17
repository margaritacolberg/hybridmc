#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# diff_s_bias.py determines the entropy difference for each transition
# whose initial and final states differ by one bond; for cases where the
# transition has a staircase potential, the entropy for each step of the
# staircase is placed into a list, from which the overall entropy of the
# staircase is found

import csv
import glob
import h5py
import numpy as np
import re
from itertools import combinations


def get_diff_sbias(out_csv='diff_s_bias.csv'):
    src = 'hybridmc_*.h5'

    bits = []
    s_bias = []
    duplicate_s_bias = {}

    for file_path in glob.glob(src):
        match = re.search(r'_(?P<bits_in>[01]+)_(?P<bits_out>[01]+)(?:_(?P<stair>[0-9]+\.[0-9]+))?\.h5$',
                          file_path)
        bits_in = match.group('bits_in')
        bits_out = match.group('bits_out')

        with h5py.File(file_path, 'r') as f:
            get_s_bias = f['s_bias'][0]
            if (bits_in, bits_out) not in bits:
                s_bias.append(get_s_bias)
            else:
                # extract entropy from output file for all steps of staircase
                # potential except last step (each step of staircase has the
                # same set of bits associated with it, hence why it is
                # 'duplicate')
                duplicate_s_bias.setdefault((bits_in, bits_out),
                                            []).append(get_s_bias)

        if (bits_in, bits_out) not in bits:
            bits.append((bits_in, bits_out))

    for i in range(len(bits)):
        if bits[i] in duplicate_s_bias.keys():
            duplicate_s_bias[bits[i]].append(s_bias[i])

    for key, value in duplicate_s_bias.items():
        exp_s = []

        for i in range(1, len(value) + 1):
            [exp_s.append(np.exp(sum(j))) for j in combinations(value, i)]

        for i in range(len(bits)):
            if key == bits[i]:
                s_bias[i] = np.log(np.sum(exp_s))

    output = []
    for i in range(len(bits)):
        output.append([bits[i][0], bits[i][1], s_bias[i]])

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)


if __name__ == '__main__':
    import os

    os.chdir('../../examples/test_with_wanglandau_also_cutoff')
    get_diff_sbias(out_csv='diff_s_bias_sort.csv')
