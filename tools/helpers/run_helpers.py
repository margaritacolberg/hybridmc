#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause


from .init_json_helpers import *
from multiprocessing import Process, Queue


def init_json(args):
    with open(args["json"], 'r') as input_json:
        data = json.load(input_json)

    # make copy of nonlocal bonds provided from json file for processing
    nonlocal_bonds = data['nonlocal_bonds']
    nbonds = len(nonlocal_bonds)
    nstates = 2 ** nbonds

    # use make rc tuple to make each nonlocal bond to be a triplet with rc included in case rc value not already
    # included as a triplet for each bond pair list
    nonlocal_bonds = make_rc_tuple(nonlocal_bonds, data["rc"])

    # sort the nonlocal_bonds list
    nonlocal_bonds = sort_triplet(nonlocal_bonds)
    nonlocal_bonds.sort()

    # set the data json files nonlocal bonds list to processes nonlocal bonds list
    data['nonlocal_bonds'] = nonlocal_bonds

    # set the transient bonds list to be the nonlocal bonds list initially
    data['transient_bonds'] = data["nonlocal_bonds"]

    in_queue = Queue()
    out_queue = Queue()

    worker = []
    layer_args = (data, in_queue, out_queue, args["seed_increment"], args["WL_sbias"])

    for _ in range(args["nproc"]):
        p = Process(target=run_layer, args=layer_args)
        p.start()
        worker.append(p)

    count = 0

    # initially no bonds
    out_queue.put(((False,) * nbonds, None))

    bits_seen = {}
    while len(bits_seen) != nstates:
        item = out_queue.get()

        bits_in, input_hdf5 = item

        # check if bond pattern already exists inside worklist regardless of
        # what hdf5 input file it is associated with; only append unique bond
        # patterns to worklist
        if bits_in in bits_seen:
            continue

        else:
            bits_seen[bits_in] = True

        for i in range(nbonds):
            # count is to initialize seed of random number generator
            count += 1

            # if two beads are bonded,
            if bits_in[i] is True:
                # skip; do not break the bond
                continue

            else:
                in_queue.put((i, bits_in, count, input_hdf5))

    # signal to run_layer that no more items are left to be added to queue
    # (native state has been reached)
    for _ in range(len(worker)):
        in_queue.put(None)

    for p in worker:
        p.join()