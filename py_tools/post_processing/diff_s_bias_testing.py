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
import asyncio
import threading
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count

if __name__ == '__main__' and (__package__ is None or __package__ == ''):
    from py_tools.helpers.bootstrap_helpers import ConfigEntropyDiffBoot, StairConfigEntropyDiffBoot
else:
    from ..helpers.bootstrap_helpers import ConfigEntropyDiffBoot, StairConfigEntropyDiffBoot

import time


def measure_execution_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Time taken for '{func.__name__}' function: {execution_time:.4f} seconds")
        return result

    return wrapper


# Example usage:

async def __async_get_diff_sbias(out_csv):
    # Initializations
    output = []  # The list to write out to csv file, each element is one row of info for each simulation
    sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

    # Prepare tasks for all simulations
    tasks = [process_simulation(sim, is_stair=False) for sim in sims]
    tasks += [process_simulation(sim, is_stair=True) for sim in stair_sims]

    # Execute tasks concurrently and collect results
    results = await asyncio.gather(*tasks)
    output.extend(results)

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

    print("Done writing diff s_bias output")


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


async def process_simulation(simulation_name, is_stair):
    thread = threading.current_thread()
    print(f'name={thread.name}, daemon={thread.daemon}')
    if is_stair:
        # For staircased simulations
        result = await asyncio.to_thread(StairConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output)

    else:
        # For normal simulations
        result = await asyncio.to_thread(ConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output)

    return result


def fake_sync_process_simulation(simulation_name, is_stair):
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        return loop.run_until_complete(process_simulation(simulation_name, is_stair=is_stair))
    finally:
        loop.close()


# Asynchronous function to process a simulation
def _sync_process_simulation(simulation_name, is_stair):
    thread = threading.current_thread()
    print(f'name={thread.name}, daemon={thread.daemon}')

    from os import getpid
    print(f"Current Process ID (PID): {getpid()}")

    if is_stair:
        # For staircased simulations
        return StairConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output()

    else:
        # For normal simulations
        return ConfigEntropyDiffBoot(simulation_name).get_diff_sbias_output()


# Synchronous wrapper for process_simulation
@measure_execution_time
def diff_s_bias_multi(out_csv):
    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        # Initializations
        output = []  # The list to write out to csv file, each element is one row of info for each simulation
        sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

        # Prepare tasks for all simulations
        tasks = [executor.submit(_sync_process_simulation, sim, is_stair=False) for sim in sims]
        tasks += [executor.submit(_sync_process_simulation, sim, is_stair=True) for sim in stair_sims]

        # Collect results from completed tasks
        results = [task.result() for task in as_completed(tasks)]
        output.extend(results)

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

        print("Done writing diff s_bias output")

@measure_execution_time
def diff_s_bias_multi(out_csv):
    with ThreadPoolExecutor(max_workers=cpu_count()) as executor:
        # Initializations
        output = []  # The list to write out to csv file, each element is one row of info for each simulation
        sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

        # Prepare tasks for all simulations
        tasks = [executor.submit(_sync_process_simulation, sim, is_stair=False) for sim in sims]
        tasks += [executor.submit(_sync_process_simulation, sim, is_stair=True) for sim in stair_sims]

        # Collect results from completed tasks
        results = [task.result() for task in as_completed(tasks)]
        output.extend(results)

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

        print("Done writing diff s_bias output")


@measure_execution_time
def get_diff_sbias_threaded_async(out_csv='diff_s_bias.csv'):
    # Create a ThreadPoolExecutor
    with ThreadPoolExecutor() as executor:
        # Initializations
        output = []  # The list to write out to csv file, each element is one row of info for each simulation
        sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

        # Prepare tasks for all simulations
        tasks = [executor.submit(fake_sync_process_simulation, sim, is_stair=False) for sim in sims]
        tasks += [executor.submit(fake_sync_process_simulation, sim, is_stair=True) for sim in stair_sims]

        # Collect results from completed tasks
        results = [task.result() for task in as_completed(tasks)]
        output.extend(results)

        with open(out_csv, 'w') as output_csv:
            writer = csv.writer(output_csv)
            output.sort()  # sort the list by the bits in column
            writer.writerows(output)

        print("Done writing diff s_bias output")


@measure_execution_time
def get_diff_sbias_async(out_csv='diff_s_bias.csv'):
    asyncio.run(__async_get_diff_sbias(out_csv))


@measure_execution_time
def diff_s_bias_sync(out_csv):
    # Initializations
    output = []  # The list to write out to csv file, each element is one row of info for each simulation
    sims, stair_sims = classify_sims('hybridmc_*.h5')  # The normal and staircase simulations

    # Prepare tasks for all simulations
    output = [_sync_process_simulation(sim, is_stair=False) for sim in sims]
    output += [_sync_process_simulation(sim, is_stair=True) for sim in stair_sims]

    with open(out_csv, 'w') as output_csv:
        writer = csv.writer(output_csv)
        output.sort()  # sort the list by the bits in column
        writer.writerows(output)

        print("Done writing diff s_bias output")


def main(out_csv='diff_s_bias_with_error.csv'):
    diff_s_bias_sync(out_csv)
    #get_diff_sbias_async(out_csv)
    # get_diff_sbias_threaded(out_csv)
    #get_diff_sbias_threaded_async(out_csv)
    #diff_s_bias_multi(out_csv)


if __name__ == '__main__':
    for _ in range(10):
        main()
