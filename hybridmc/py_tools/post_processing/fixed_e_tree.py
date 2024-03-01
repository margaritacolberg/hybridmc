#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# fixed_e_tree.py packages configurations, entropies, and energies of all
# states of the protein into .dat and .mtx files for creating tree diagrams
# from fixed_e data
#
# example of how to run:
# python ../tools/fixed_e_tree.py crambin_fixed_e

import matrix_element

import argparse
import csv
import glob


def main(args):
    src = 'hybridmc_*.csv'

    csv_in = 'avg_s_bias.csv'
    s_bias = matrix_element.get_s_bias(csv_in)

    config = []
    output_1 = []
    output_2 = []
    for file_path in glob.glob(src):
        bits_i, bits_j, t_ind, int_i, int_j = \
                matrix_element.get_state_data(file_path)

        e_in = -bin(int_i).count('1')

        if int_i not in config:
            output_1.append([int_i, e_in, s_bias[int_i], bits_i])
            config.append(int_i)

        nbits = len(bits_i)
        if int_j == (2**nbits - 1) and int_j not in config:
            e_out = -bin(int_j).count('1')
            output_1.append([int_j, e_out, s_bias[int_j], bits_j])
            config.append(int_j)

        dS_of_bond_on = s_bias[int_j] - s_bias[int_i]
        dS_of_bond_off = s_bias[int_i] - s_bias[int_j]

        inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)

        output_2.append([int_i, int_j, outer_fpt, dS_of_bond_on])
        output_2.append([int_j, int_i, inner_fpt, dS_of_bond_off])

    output_1 = sorted(output_1, key=lambda x: x[0])
    output_2 = sorted(output_2, key=lambda x: (x[0], x[1]))

    dat_file = '{}.dat'.format(args.name)
    with open(dat_file, 'w') as output_csv:
        writer = csv.writer(output_csv, delimiter=' ')
        writer.writerows(output_1)

    mtx_file = '{}.mtx'.format(args.name)
    with open(mtx_file, 'w') as output_csv:
        writer = csv.writer(output_csv, delimiter=' ')
        writer.writerows(output_2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('name', help='name of output files')

    args = parser.parse_args()

    main(args)
