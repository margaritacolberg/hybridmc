#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import copy
import json
import os
import subprocess
import sys
from multiprocessing import Process, Queue


def main(args):
    with open(args.json, 'r') as input_json:
        data = json.load(input_json)

    nonlocal_bonds = data['nonlocal_bonds']
    nbonds = len(nonlocal_bonds)
    nstates = 2**nbonds

    nonlocal_bonds = [sorted(el) for el in nonlocal_bonds]
    nonlocal_bonds.sort()

    # use make rc tuple to make each nonlocal bond to be a triplet with rc included
    data['nonlocal_bonds'] = make_rc_tuple(nonlocal_bonds, data["rc"])
    data['transient_bonds'] = data["transient_bonds"]

    in_queue = Queue()
    out_queue = Queue()

    worker = []
    layer_args = (nonlocal_bonds, data, in_queue, out_queue, args.exe,
            args.seed_increment, args.rc, args.bp)

    for i in range(args.nproc):
        p = Process(target=run_layer, args=layer_args)
        p.start()
        worker.append(p)

    count = 0

    # initially no bonds
    out_queue.put(((False,)*nbonds, None))

    bits_seen = {}
    while len(bits_seen) != nstates:
        item = out_queue.get()

        bits_in, input_hdf5 = item

        # check if bond pattern already exists inside worklist regardless of
        # what hdf5 input file it is associated with; only append unique bond
        # patterns to worklist
        if bits_in in bits_seen:
            continue

        bits_seen[bits_in] = True

        for i in range(nbonds):
            # count is to initialize seed of random number generator
            count += 1

            # if two beads are bonded,
            if bits_in[i]:
                # skip; do not break the bond
                continue

            in_queue.put((i, bits_in, count, input_hdf5))

    # signal to run_layer that no more items are left to be added to queue
    # (native state has been reached)
    for i in range(len(worker)):
        in_queue.put(None)

    for p in worker:
        p.join()


def run_layer(nonlocal_bonds, common_data, in_queue, out_queue, exe,
        seed_increment, rc, bp):
    while True:
        item = in_queue.get()
        if item is None:
            break

        i, bits_in, count, input_hdf5 = item

        # copy of initial state for every nonbonded pair of beads
        bits_out = list(bits_in)
        # flip bit to form bond
        bits_out[i] = True
        bits_out = tuple(bits_out)

        bonds_in = bits_to_bonds(bits_in, nonlocal_bonds)
        layer = len(bonds_in)

        # optional staircase potential
        if bp is not None and nonlocal_bonds[i] == bp:
            for j in range(len(rc)):
                data = copy.deepcopy(common_data)

                # if there is at least one permanent bond,
                if layer > 0:
                    data['p_rc'] = common_data['rc']

                if j > 0:
                    data['stair_bonds'] = [nonlocal_bonds[i]]
                    data['stair'] = rc[j-1]

                if j < len(rc)-1:
                    output_name = 'hybridmc_{}_{}_{}_{}'.format(layer,
                            format_bits(bits_in), format_bits(bits_out), rc[j])
                    hdf5_name = '{}.h5'.format(output_name)
                else:
                    output_name = 'hybridmc_{}_{}_{}'.format(layer,
                            format_bits(bits_in), format_bits(bits_out))
                    hdf5_name = '{}.h5'.format(output_name)

                data['rc'] = rc[j]
                run_sim(nonlocal_bonds[i], data, exe, seed_increment,
                        input_hdf5, hdf5_name, output_name, bits_in, bits_out,
                        bonds_in, count)
                input_hdf5 = hdf5_name
        else:
            output_name = 'hybridmc_{}_{}_{}'.format(layer,
                    format_bits(bits_in), format_bits(bits_out))
            hdf5_name = '{}.h5'.format(output_name)
            run_sim(nonlocal_bonds[i], common_data, exe, seed_increment,
                    input_hdf5, hdf5_name, output_name, bits_in, bits_out,
                    bonds_in, count)

        out_queue.put((bits_out, hdf5_name))


def run_sim(nonlocal_bonds_i, data, exe, seed_increment, input_hdf5, hdf5_name,
        output_name, bits_in, bits_out, bonds_in, count):
    data['transient_bonds'] = [nonlocal_bonds_i]
    data['permanent_bonds'] = bonds_in
    data['config_in'] = int(format_bits(bits_in), 2)
    data['config_out'] = int(format_bits(bits_out), 2)
    data['seeds'] = [count, seed_increment]

    if os.path.exists(hdf5_name):
        return

    json_name = '{}.json'.format(output_name)
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # for layer = 1 or greater,
    command = [exe, json_name, hdf5_name]
    if not input_hdf5 is None:
        command += ['--input-file', input_hdf5]

    print(command)
    sys.stdout.flush()

    log_name = '{}.log'.format(output_name)
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                stderr=subprocess.STDOUT)


def bits_to_bonds(bits, nonlocal_bonds):
    bits = [i for i, bit in enumerate(bits) if bit]
    bonds = [nonlocal_bonds[i] for i in bits]
    return bonds


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))


def make_rc_tuple(nonlocal_bonds, rc):
    """
    Function to produce nonlocal bonds list with each nonlocal bond list element containing the rc as the third
    element

    args:
    nonlocal_bonds -- List of Lists containing the indices of beads bonded to each other, for example [ [1, 9], [4, 15]].
    These elements may potentially have the bond length (rc) as the third element as well or not,
    for example [[1, 9, 2.3], [2, 6]].

    rc -- Default rc to be appended to each nonlocal_bonds element in case they do not have this information

    returns:
    The same list but with each list element having the bond length (rc) as their third element.

    """

    # loop through each bonded bead pair
    for el in nonlocal_bonds:

        # check if rc for the bond is already given
        if len(el) != 3:

            # if not then append rc
            el.append(rc)

    return nonlocal_bonds


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='master json input file')
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('seed_increment', type=int,
            help='amount to increment seed with every run')
    parser.add_argument('--nproc', type=int, default=os.cpu_count(),
            help='number of worker processes')
    parser.add_argument('--rc', nargs='+', type=float,
            help='list of rc if potential is step')
    parser.add_argument('--bp', nargs=2, type=int,
            help='bead pair whose potential is step')

    args = parser.parse_args()

    main(args)


