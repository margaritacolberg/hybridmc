#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# diff_s_bias.py determines the entropy difference for each transition
# whose initial and final states differ by one bond; for cases where the
# transition has a staircase potential, the entropy for each step of the
# staircase is placed into a list, from which the overall entropy of the
# staircase is found

import csv
import glob
import re

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from py_tools.helpers.bootstrap_helpers import ConfigBoot, StairConfigBoot
else:
    from ..helpers.bootstrap_helpers import ConfigBoot, StairConfigBoot


def get_diff_sbias_pair(base_file):
    src = f'{base_file}*.h5'
    bits = []
    s_bias = []
    duplicate_s_bias = {}
    for file_path in glob.glob(src):
        match = re.search(r'_(?P<bits_in>[01]+)_(?P<bits_out>[01]+)(?:_(?P<stair>[0-9]+\.[0-9]))',
                          file_path)
        bits_in = match.group('bits_in')
        bits_out = match.group('bits_out')


def get_diff_sbias(out_csv='diff_s_bias.csv'):
    src = 'hybridmc_*.h5'

    # Initializations
    output = []  # The list to write out to csv file, each element is one row of info for each simulation
    sims = set()  # The set of normal simualations
    stair_sims = set()  # set of staircase simualtions

    for file_path in glob.glob(src):
        match = re.search(r'_(?P<bits_in>[01]+)_(?P<bits_out>[01]+)(?:_(?P<stair>[0-9]+\.[0-9]+))?\.h5$',
                          file_path)

        # check if this is a staircase intermediate simualation
        stair = match.group('stair')
        simulation_name = ConfigBoot.get_base_simulation_id(file_path)

        # add the simulation name to the stair_sims set and the normal ones to the sims set
        if stair:
            stair_sims.add(simulation_name)
        else:
            sims.add(simulation_name)

        # remove all the stair sim names from the sim names set
        #sims = sims - stair_sims

        # in case this simulation name is a staircased one use StairConfig routine
        for simulation_name in stair_sims:
            get_config_info = StairConfigBoot(simulation_name=simulation_name)
            output.append(get_config_info.get_diff_sbias_output())

        # In case it is a normal simulation, use Config routine
        for simulation_name in sims:
            get_config_info = ConfigBoot(simulation_name=simulation_name)
            output.append(get_config_info.get_diff_sbias_output())

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

    print("Done writing diff s_bias output")


if __name__ == '__main__':
    get_diff_sbias(out_csv='diff_s_bias_with_error.csv')
