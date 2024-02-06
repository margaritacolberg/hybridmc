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
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from py_tools.helpers.bootstrap_helpers import ConfigEntropyDiffBoot, StairConfigEntropyDiffBoot
else:
    from ..helpers.bootstrap_helpers import ConfigEntropyDiffBoot, StairConfigEntropyDiffBoot


# Example usage:


def classify_sims(src):
    sims = set()  # The set of normal simualations
    stair_sims = set()  # set of staircase simualtions
    for file_path in glob.glob(src):
        match = re.search(r'(?P<simulation_name>hybridmc_(\d+)_([01]+)_([01]+))(?:_(?P<stair>[0-9]+\.[0-9]+))?\.h5$',
                          file_path)

        # add the simulation name to the stair_sims set and the normal ones to the sims set if it is a staircase h5 file
        if match.group('stair'):
            stair_sims.add(match.group('simulation_name'))
        else:
            sims.add(match.group('simulation_name'))

    # remove all the stair sim names from the sim names set
    sims = sims - stair_sims
    return sims, stair_sims


# Synchronous wrapper for process_simulation
def process_simulation(simulation_name, is_stair):
    print(f"Processing {simulation_name} ...")
    return StairConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output() if is_stair else \
        ConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output()


def multi_diff_s_bias(out_csv):
    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:

        # Initializations
        output = []  # The list to write out to csv file, each element is one row of info for each simulation
        sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

        # Prepare tasks for all simulations
        tasks = [executor.submit(process_simulation, sim, is_stair=False) for sim in sims]
        tasks += [executor.submit(process_simulation, sim, is_stair=True) for sim in stair_sims]

        # Collect results from completed tasks
        results = [task.result() for task in as_completed(tasks)]
        output.extend(results)

    # write out the results to csv file
    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

    print("Done writing diff s_bias output")


def get_diff_sbias(out_csv='diff_s_bias_with_error.csv'):
    multi_diff_s_bias(out_csv)


if __name__ == '__main__':
    get_diff_sbias()

