#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# testing.py is used for debugging a single transition for which the simulation
# compiles successfully but crashes when it starts running
#
# example of how to run:
# python ../tools/testing.py hybridmc_1_100_110.json ../release/hybridmc --hdf5 hybridmc_0_000_100.h5

import argparse
import subprocess
import sys


def main(args):
    hdf5_name = 'test.h5'

    command = ['time', args.exe, args.json, hdf5_name]
    if not args.hdf5 is None:
        command += ['--input-file', args.hdf5]

    print(command)
    sys.stdout.flush()

    log_name = 'test.log'
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                stderr=subprocess.STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('json', help='json input file')
    parser.add_argument('exe', help='hybridmc executable')
    parser.add_argument('--hdf5', help='hdf5 input file')

    args = parser.parse_args()

    main(args)
