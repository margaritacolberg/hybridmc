#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import matrix_element

import csv
import os
import re


def main():
    output = []
    for file_path in os.listdir('.'):
        if re.search('hybridmc_[0-9]+_([01]+)_([01]+).csv', file_path):
            bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(file_path)
            inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)
            layer = str(bits_i).count('1')
            output.append([bits_i, bits_j, layer, inner_fpt, outer_fpt])

    output.sort(key=lambda x: x[2])
    csv_out = 'master_fpt.csv'

    output.insert(0, ['state i bits', 'state j bits', 'layer', 'inner fpt',
        'outer fpt'])

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


if __name__ == '__main__':
    main()
