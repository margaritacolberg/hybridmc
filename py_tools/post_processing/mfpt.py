#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# mfpt.py determines the mean first passage times using spline interpolation
# and their errors using the bootstrap method; the outer first passage time
# (fpt) corresponds to the first time at which two beads at a distance r > rc
# diffuse to rc, and the inner fpt corresponds to beads located at r < rc
#
#
# note that in the above example, hybridmc_0_000_001.json is the json file for
# the last step of the staircase (ie. where the width of the step is between rh
# and rc)

import csv
import os
import re
import glob
from . import matrix_element
import multiprocessing

from ..helpers.mfpt_helpers import fpt_write


def get_mfpt_serial(rewrite=False):

    # iterate over all json files in the current directory
    src = '*.json'
    names = []
    for file_path in glob.glob(src):
        file_name = os.path.basename(file_path)
        name = os.path.splitext(file_name)[0]

        if os.path.exists(f"{name}.csv") and not rewrite:
            continue

        names.append(name)

    # write the mfpt for each json file found serially

    for name in names:
         fpt_write(name)


def get_mfpt(rewrite=False):

    # iterate over all json files in the current directory
    src = '*.json'
    names = []
    for file_path in glob.glob(src):
        file_name = os.path.basename(file_path)
        name = os.path.splitext(file_name)[0]

        if os.path.exists(f"{name}.csv") and not rewrite:
            continue

        names.append(name)

    # write the mfpt for each json file found in parallel
    with multiprocessing.Pool() as pool:
       pool.map(fpt_write, names)


def compile_mfpts():
    """
    Function to compile all the mfpt csv files in the current directory into one csv file

    Returns
    -------
    None, but writes a csv file with the compiled mfpt data

    """
    output = []
    files = os.listdir()
    for file_path in files:

        if re.search('hybridmc_[0-9]+_([01]+)_([01]+).csv', file_path):

            bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(file_path)
            inner_fpt, outer_fpt, inner_std, outer_std = matrix_element.get_fpt(file_path, t_ind)
            layer = file_path.split('_')[1]
            output.append([bits_i, bits_j, layer, inner_fpt, outer_fpt, inner_std, outer_std])

    output.sort(key=lambda x: (x[2],int(x[0],2),int(x[1],2)))
    csv_out = 'mfpt.csv'

    output.insert(0, ['state i bits', ' state j bits', ' layer', ' inner fpt',
                      ' outer fpt', ' inner std', ' outer std'])

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)
