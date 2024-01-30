#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# min_e.py determines the average time it takes for a protein to fold to the
# native state, starting from the unfolded state, where the probability of
# encountering the native state is kept high relative to the unfolded state
# (the minimum energy of each intermediate state is calculated)
#
# example of how to run:
# python ../tools/min_e.py ../examples/crambin.json -60.0 min_e.csv

import matrix_element
import minimize
import argparse
import csv
import glob
import json
import nlopt
import numpy as np


def main(args):
    src = 'hybridmc_*.csv'

    with open(args.json, 'r') as input_json:
        json_data = json.load(input_json)

    nbonds = json_data['max_nbonds']
    nstates = 2**nbonds
    ne = nstates-2
    beta = 1
    e_unfolded = 0.0

    csv_in = 'avg_s_bias.csv'
    s_bias = matrix_element.get_s_bias(csv_in)
    s_bias = np.array(s_bias)

    e = np.zeros(ne)
    for i in range(1, nstates-2):
        bits_i = '{0:b}'.format(i)
        e[i] = -bits_i.count('1')

    ub = np.zeros(ne)
    lb = np.zeros(ne)

    prob_ratio = 0.01
    lb[:nstates-2] = s_bias[1:nstates-1] + args.e_native - np.log(prob_ratio)

    df_wrapped = lambda x0: df(src, x0, beta, e_unfolded, args.e_native,
            s_bias, nstates)
    f_wrapped = lambda x0: f(src, x0, beta, e_unfolded, args.e_native, s_bias,
            nstates)

    # note that the last function evaluation during minimization does not
    # necessarily produce the final optimized values in e
    e = minimize.minimize(f_wrapped, e, method=nlopt.GN_DIRECT_L,
            jac=df_wrapped, ub=ub, lb=lb)
    e = minimize.minimize(f_wrapped, e, method=nlopt.LD_LBFGS, jac=df_wrapped,
            ub=ub, lb=lb)

    P_0 = np.exp(s_bias[0])
    P_native = np.exp(-beta*args.e_native + s_bias[nstates-1])
    P_norm = P_0 + np.sum(np.exp(-beta*e[:] + s_bias[1:nstates-1])) + P_native

    P = []
    P.append(P_0 / P_norm)
    P[1:nstates-1] = np.exp(-beta*e[:] + s_bias[1:nstates-1]) / P_norm
    P.append(P_native / P_norm)

    output = []

    min_e = f(src, e, beta, e_unfolded, args.e_native, s_bias, nstates)
    output.append([min_e])

    for i in range(nstates):
        if i == 0:
            print('state i = {} has E = {}, S = {}, and P = {}'.format(i, -0.0,
                s_bias[i], P[i]))
            output.append([i, -0.0, s_bias[i], P[i]])
        elif i == nstates-1:
            print('state i = {} has E = {}, S = {}, and P = {}'.format(i,
                args.e_native, s_bias[i], P[i]))
            output.append([i, args.e_native, s_bias[i], P[i]])
        else:
            print('state i = {} has E = {}, S = {}, and P = {}'.format(i,
                e[i-1], s_bias[i], P[i]))
            output.append([i, e[i-1], s_bias[i], P[i]])

    with open(args.csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def f(src, e, beta, e_unfolded, e_native, s_bias, nstates):
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

        if int_i == 0:
            K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                    e_unfolded, e[int_j-1], s_bias[int_i], s_bias[int_j],
                    t_ind)
        elif int_j == nstates-1:
            K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                    e[int_i-1], e_native, s_bias[int_i], s_bias[int_j],
                    t_ind)
        else:
            K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                    e[int_i-1], e[int_j-1], s_bias[int_i], s_bias[int_j],
                    t_ind)

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


def df(src, e, beta, e_unfolded, e_native, s_bias, nstates):
    dmfpt = []

    for k in range(1, nstates-1):
        state_i = []
        state_j = []
        K_ji = []
        K_ij = []
        dK_ji = []
        dK_ij = []

        for file_path in glob.glob(src):
            bits_i, bits_j, t_ind, int_i, int_j = \
                    matrix_element.get_state_data(file_path)

            state_i.append(int_i)
            state_j.append(int_j)

            if int_i == 0:
                K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                        e_unfolded, e[int_j-1], s_bias[int_i], s_bias[int_j],
                        t_ind)
                dK_elem_ji, dK_elem_ij = matrix_element.dK_elem(file_path,
                        beta, e_unfolded, e[int_j-1], s_bias[int_i],
                        s_bias[int_j], int_i, int_j, k, t_ind)
            elif int_j == nstates-1:
                K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                        e[int_i-1], e_native, s_bias[int_i], s_bias[int_j],
                        t_ind)
                dK_elem_ji, dK_elem_ij = matrix_element.dK_elem(file_path,
                        beta, e[int_i-1], e_native, s_bias[int_i],
                        s_bias[int_j], int_i, int_j, k, t_ind)
            else:
                K_elem_ji, K_elem_ij = matrix_element.K_elem(file_path, beta,
                        e[int_i-1], e[int_j-1], s_bias[int_i], s_bias[int_j],
                        t_ind)
                dK_elem_ji, dK_elem_ij = matrix_element.dK_elem(file_path,
                        beta, e[int_i-1], e[int_j-1], s_bias[int_i],
                        s_bias[int_j], int_i, int_j, k, t_ind)

            K_ji.append(K_elem_ji)
            K_ij.append(K_elem_ij)
            dK_ji.append(dK_elem_ji)
            dK_ij.append(dK_elem_ij)

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

        dK = np.zeros((nstates, nstates))

        for i in range(len(state_i)):
            for j in range(len(state_i)):
                if i == state_i[j]:
                    dK[state_j[j]][state_i[j]] = dK_ji[j]
                    dK[state_i[j]][state_j[j]] = dK_ij[j]

        for i in range(nstates):
            dK[i][i] = -np.sum(dK[:, i])

        dK_tilde = np.delete(dK, nstates-1, axis=0)
        dK_tilde = np.delete(dK_tilde, nstates-1, axis=1)

        matrix_1 = np.dot(dK_tilde, K_tilde_inv)
        matrix_2 = np.dot(-K_tilde_inv, matrix_1)

        dmfpt.append(-np.sum(matrix_2[:, 0]))

    return np.array(dmfpt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('e_native', type=float, help='energy of native state')
    parser.add_argument('csv', help='csv output file')

    args = parser.parse_args()

    main(args)
