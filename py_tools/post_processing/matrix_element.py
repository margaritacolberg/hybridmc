#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# matrix_element.py calcuates each element of the transition rate matrix and
# its derivative with respect to energy

import csv
import math
import re


def K_elem(csv_in, beta, e_i, e_j, s_i, s_j, t_ind):
    inner_fpt, outer_fpt, inner_std, outer_std = get_fpt(csv_in, t_ind)

    # the final state j has one more bond than the initial state i
    s_bias_diff = s_i - s_j
    eps = e_i - e_j

    P_i_div_P_j = math.exp((-beta*eps) + s_bias_diff)
    K_ji_inv = P_i_div_P_j*inner_fpt + outer_fpt
    K_ji = 1/K_ji_inv
    K_ij = P_i_div_P_j*K_ji

    return K_ji, K_ij


def dK_elem(csv_in, beta, e_i, e_j, s_i, s_j, i, j, k, t_ind):
    inner_fpt, outer_fpt, inner_std, outer_std = get_fpt(csv_in, t_ind)

    s_bias_diff = s_i - s_j
    eps = e_i - e_j

    P_i_div_P_j = math.exp((-beta*eps) + s_bias_diff)
    K_ji_inv = P_i_div_P_j*inner_fpt + outer_fpt
    K_ji = 1/K_ji_inv
    K_ij = P_i_div_P_j*K_ji

    if k == j and k != i:
        # dK_ji / de_j
        dK_ji = -beta*P_i_div_P_j*inner_fpt*K_ji**2
        # dK_ij / de_j
        dK_ij = (beta*K_ij) + (P_i_div_P_j*dK_ji)
    elif k != j and k == i:
        # dK_ji / de_i
        dK_ji = beta*P_i_div_P_j*inner_fpt*K_ji**2
        # dK_ij / de_i
        dK_ij = -(beta*K_ij) + (P_i_div_P_j*dK_ji)
    else:
        dK_ji = 0.0
        dK_ij = 0.0

    return dK_ji, dK_ij


def get_fpt(csv_in, t_ind):
    with open(csv_in, 'r') as input_csv:
        data = csv.reader(input_csv, delimiter=',')

        for row in data:
            if int(row[0]) == t_ind:
                inner_fpt = float(row[1])
                outer_fpt = float(row[2])
                inner_std = float(row[3])
                outer_std = float(row[4])

    return inner_fpt, outer_fpt, inner_std, outer_std


def get_s_bias(csv_in):
    s_bias = []
    with open(csv_in, 'r') as input_csv:
        csv_data = csv.reader(input_csv, delimiter=',')

        for row in csv_data:
            # skip first column in avg_s_bias.csv, which is bit pattern of
            # configuration
            s_bias.append(float(row[1]))

    return s_bias


def get_state_data(file_path):
    print('input csv:', file_path)

    # get bit pattern of states i and j from title of each csv
    bits_from_csv = re.search(r'_([01]+)_([01]+)\.csv$', file_path)

    if bits_from_csv is None:
        return

    bits_i, bits_j = bits_from_csv.groups()

    # bond that differs between bit patterns of i and j is the transient
    # bond; get its index
    t_ind = [i for i in range(len(bits_i)) if bits_i[i] != bits_j[i]][0]

    int_i = int(bits_i, base=2)
    int_j = int(bits_j, base=2)

    return bits_i, bits_j, t_ind, int_i, int_j
