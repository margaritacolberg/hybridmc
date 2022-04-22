#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# diff_s_bias.py determines the entropy difference for each transition whose
# initial and final states differ by one bond
#
# example of how to run:
# python ../tools/diff_s_bias.py ../examples/crambin.json

import argparse
import csv
import h5py
import json


def main(args):
    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    nonlocal_bonds = data['nonlocal_bonds']
    nbonds = len(nonlocal_bonds)

    # initially no bonds
    work_list = [[False]*nbonds]

    output = []

    while work_list:
        bits_in = work_list.pop(0)

        for i in range(nbonds):
            # if two beads are bonded,
            if bits_in[i]:
                # skip; do not break the bond
                continue

            # copy of initial state for every nonbonded pair of beads
            bits_out = bits_in.copy()
            # flip bit to form bond
            bits_out[i] = True

            bonds_in = bits_to_bonds(bits_in, nonlocal_bonds)
            layer = len(bonds_in)

            output_name = 'hybridmc_{}_{}_{}'.format(layer,
                    format_bits(bits_in), format_bits(bits_out))
            hdf5_name = '{}.h5'.format(output_name)

            config_in = format_bits(bits_in)
            config_out = format_bits(bits_out)

            with h5py.File(hdf5_name, 'r') as f:
                output.append((config_in, config_out, float(f['s_bias'][0])))

            # check if bond pattern already exists inside worklist;
            # only append unique bond patterns to worklist
            if not bits_out in work_list:
                work_list.append(bits_out)

    with open('diff_s_bias.csv', 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def bits_to_bonds(bits, nonlocal_bonds):
    bits = [i for i, bit in enumerate(bits) if bit]
    bonds = [nonlocal_bonds[i] for i in bits]
    return bonds


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')

    args = parser.parse_args()

    main(args)
