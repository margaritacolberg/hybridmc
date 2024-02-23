#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# fixed_e.py determines the average time it takes for a protein to fold to the
# native state, starting from the unfolded state, where all nonlocal bonds are
# kept at the same energy while varying beta*
#
# example of how to run:
# python ../tools/fixed_e.py ../examples/crambin.json 5.5 fixed_e.csv

import matrix_element
import argparse
import csv
import glob
import json
import numpy as np


def main(args):
    src = 'hybridmc_*.csv'

    with open(args.json, 'r') as input_json:
        json_data = json.load(input_json)

    nbonds = json_data['max_nbonds']
    nstates = 2**nbonds

    csv_in = 'avg_s_bias.csv'
    s_bias = matrix_element.get_s_bias(csv_in)
    s_bias = np.array(s_bias)

    P = []
    for i in range(nstates):
        bits_i = '{0:b}'.format(i)
        e_i = -bits_i.count('1')
        P.append(np.exp(-args.beta*e_i + s_bias[i]))

    P_norm = np.sum(P)
    P = P / P_norm

    output = []

    fixed_e = f(src, args.beta, s_bias, nstates)
    output.append([fixed_e])

    for i in range(nstates):
        output.append([i, args.beta, s_bias[i], P[i]])

    with open(args.csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def f(src, beta, s_bias, nstates):
    # the final state j has one more bond than the initial state i
    state_i = []
    state_j = []
    K_ji = []
    K_ij = []

    for file_path in glob.glob(src):
        bits_i, bits_j, t_ind, int_i, int_j = \
                matrix_element.get_state_data(file_path)

        state_i.append(int_i)
        state_j.append(int_j)

        e_i = -bits_i.count('1')
        e_j = -bits_j.count('1')

        K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta, e_i, e_j,
                s_bias[int_i], s_bias[int_j], t_ind)

        K_ji.append(K_elem_ji)
        K_ij.append(K_elem_ij)

    K = np.zeros((nstates, nstates))

    for i in range(len(state_i)):
        for j in range(len(state_i)):
            if i == state_i[j]:
                K[state_j[j]][state_i[j]] = K_ji[j]
                K[state_i[j]][state_j[j]] = K_ij[j]

    for i in range(nstates):
        K[i][i] = -np.sum(K[:, i])

    K_tilde = np.delete(K, nstates-1, axis=0)
    K_tilde = np.delete(K_tilde, nstates-1, axis=1)
    K_tilde_inv = np.linalg.inv(K_tilde)

    # sum of first column of K_tilde_inv is the time it takes to go from state
    # 0 to the native state (sum of entries in second column is the time it
    # takes to go from state 1 to the native state, third column is from state
    # 2 to the native state, etc.)
    mfpt = -np.sum(K_tilde_inv[:, 0])

    return mfpt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('beta', type=float)
    parser.add_argument('csv', help='csv output file')

    args = parser.parse_args()

    main(args)
