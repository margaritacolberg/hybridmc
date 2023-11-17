#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# mfpt.py determines the mean first passage times using spline interpolation
# and their errors using the bootstrap method; the outer first passage time
# (fpt) corresponds to the first time at which two beads at a distance r > rc
# diffuse to rc, and the inner fpt corresponds to beads located at r < rc
#
# example of how to run:
# python ../tools/mfpt.py hybridmc_0_000_001.json hybridmc_0_000_001.h5 10 hybridmc_0_000_001.csv --layers=True
#
# note that in the above example, hybridmc_0_000_001.json is the json file for
# the last step of the staircase (ie. where the width of the step is between rh
# and rc)

import csv
import os
import re
import glob
import matrix_element
import multiprocessing
from ..helpers.mfpt_helpers import fpt_write


def process_mfpts():
    output = []
    for file_path in os.listdir('.'):
        if re.search('hybridmc_[0-9]+_([01]+)_([01]+).csv', file_path):
            bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(file_path)
            inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)
            layer = str(bits_i).count('1')
            output.append([bits_i, bits_j, layer, inner_fpt, outer_fpt])

    output.sort(key=lambda x: x[2])
    csv_out = 'mfpt.csv'

    output.insert(0, ['state i bits', 'state j bits', 'layer', 'inner fpt',
                      'outer fpt'])

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def get_mfpt():
    src = '*.json'
    names = []
    for file_path in glob.glob(src):
        file_name = os.path.basename(file_path)
        name = os.path.splitext(file_name)[0]

        json_in = name + '.json'
        hdf5_in = name + '.h5'
        csv_out = name + '.csv'

        if os.path.exists(csv_out):
            continue

        names.append((json_in, hdf5_in, 0, csv_out, False))

    with multiprocessing.Pool() as pool:
        pool.starmap(fpt_write, names)