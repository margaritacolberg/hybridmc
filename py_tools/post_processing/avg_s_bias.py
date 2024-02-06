#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# avg_s_bias.py adds the entropy differences of all transitions from the native
# state to the current state, and this sum becomes the entropy of the current
# state (the entropy of the native state is 0); whenever two (or more)
# transitions lead to the same current state, an average of the entropy values
# is taken
#
# example of how to run:
# python ../tools/avg_s_bias.py ../examples/crambin.json diff_s_bias.csv

import asyncio
import csv
from statsmodels.stats.weightstats import DescrStatsW
import numpy as np


async def find_paths(start_bitstring, termination_bitstring, current_path, all_paths):
    if start_bitstring == termination_bitstring:
        all_paths.append(current_path.copy())
        return

    for i in range(len(start_bitstring)):
        if start_bitstring[i] != termination_bitstring[i]:
            next_bitstring = start_bitstring[:i] + \
                             str(int(not int(start_bitstring[i]))) + \
                             start_bitstring[i + 1:]
            current_path.append(next_bitstring)
            await find_paths(next_bitstring, termination_bitstring, current_path, all_paths)
            current_path.pop()


async def process_path(path, transition_info):
    summed_entropy = 0
    summed_var = 0
    for i in range(len(path) - 1):
        bitstring_in = path[i]
        bitstring_out = path[i + 1]
        for transition in transition_info:
            if transition[0] == bitstring_in and transition[1] == bitstring_out:
                summed_entropy += transition[2]
                summed_var += transition[3] ** 2
                break
    return summed_entropy, summed_var


async def process_all_paths(all_paths, transition_info):
    master_entropy_var_list = []
    tasks = []
    for path in all_paths:
        task = asyncio.create_task(process_path(path, transition_info))
        tasks.append(task)
    results = await asyncio.gather(*tasks)
    for result in results:
        master_entropy_var_list.append(result)

    # Get the transpose of the entropy and weights list to convert tuples to separate arrays of the
    # entropy and its variance
    entropy_weights_data = np.array(master_entropy_var_list).T

    # calculate the weights as the inverse variance (give preference to more accurate values which have smaller
    # variance)
    entropy_weights_data[1] = 1 / entropy_weights_data[1]

    return entropy_weights_data


async def avg_s_bias_weighted(diff_s_bias_csv):
    transition_info = []  # Load transition info from CSV file

    # Simulated transition info, replace this with actual loading from CSV
    with open(diff_s_bias_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            transition_info.append((row[0], row[1], float(row[2]), float(row[3])))

    start_bitstrings = set(list(zip(*transition_info))[0])
    termination_bitstring = '1' * len(transition_info[0][0])

    output = [(termination_bitstring, 0.0, 0.0)]
    for start_bitstring in start_bitstrings:
        all_paths = []
        await find_paths(start_bitstring, termination_bitstring, [start_bitstring], all_paths)

        entropy_weights_data = await process_all_paths(all_paths, transition_info)

        d1 = DescrStatsW(entropy_weights_data[0], weights=entropy_weights_data[1])

        # find the standard error of the weighted average by taking 1 / sqrt(the sum of the weights for each path)
        std_error = 1 / np.sqrt(d1.sum_weights) if d1.sum_weights else 0

        output.append((start_bitstring, d1.mean, std_error))
        print("FINISHED " + start_bitstring)

    return sorted(output, key=lambda t: (t[0].count('1'), int(t[0], 2)))


def get_avg_sbias(diff_sbias_csv='diff_s_bias_with_error.csv', output_csv='avg_s_bias_weighted.csv'):
    entropy_error_data = asyncio.run(avg_s_bias_weighted(diff_sbias_csv))

    with open(output_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(entropy_error_data)

    print("Done writing avg s_bias output")


if __name__ == '__main__':
    get_avg_sbias()
