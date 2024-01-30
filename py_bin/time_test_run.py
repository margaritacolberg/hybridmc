from functools import wraps
import time
import json
import csv
import argparse
import shutil
import argparse
import os
from py_tools.helpers.run_helpers import init_json
from py_tools.post_processing import diff_s_bias, avg_s_bias, mfpt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return total_time

    return timeit_wrapper

def get_json_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def update_json_data(json_data, new_values):
    for key, value in new_values.items():
        json_data[key] = value
    return json_data


def save_json_data(json_data, file_path):
    with open(file_path, 'w') as f:
        json.dump(json_data, f, indent=4)


def test_runtimes(file_path, new_values):
    data = get_json_data(file_path)
    data = update_json_data(data, new_values)
    save_json_data(data, 'runtime_test.json')

    args = {}
    args['json'] = 'runtime_test.json'
    args['exe'] = "../release/hybridmc"

    return hmc(args)


@timeit
def hmc(args):
    """
    Simple function that runs main
    """

    file_name = os.path.basename(args["json"])
    dir_name = os.path.splitext(file_name)[0]

    # Change directory name to suit the new temp working directory.
    # add ../ to the path to indicate its use from a directory one more level down.
    (args['json'], args['exe']) = ("../" + path for path in (args["json"], args["exe"]))

    # Create dictionary that will have arguments passed to init_json
    init_json_args = {"json": args["json"], "seed_increment": 1, "exe": args["exe"]}

    nproc = os.cpu_count()
    if os.getenv('SLURM_CPUS_PER_TASK'):
        nproc = os.getenv('SLURM_CPUS_PER_TASK')

    init_json_args["nproc"] = nproc

    # create a temporary directory to run the simulations
    tmp_dir_name = f'{dir_name}.tmp'
    if not os.path.isdir(tmp_dir_name):
        os.mkdir(tmp_dir_name)

    # move into the temporary directory
    os.chdir(tmp_dir_name)

    # run the simulations for the layers
    init_json(init_json_args)

    # Move up from the directory with simulation results
    os.chdir("../")

    # Rename the directory -- remove the .tmp tag to show that this simulation has run completely with success
    shutil.rmtree(tmp_dir_name, ignore_errors=True, onerror=None)


def time_test(name):

    nbeads = [30, 40]
    del_t = [1, 2, 3, 4, 5]
    times = [[], [], []]
    for nbead in nbeads:
        for dt in del_t:
            run_time = test_runtimes(f'{name}.json', {'nbeads': nbead, 'del_t': dt})
            times[0].append(nbead), times[1].append(dt), times[2].append(round(run_time, 2))

    with open(f'{name}_time_test.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['nbeads', 'del_t', 'time'])
        writer.writerows(zip(*times))

    df = pd.read_csv(f'{name}_time_test.csv')
    fig, ax = plt.subplots()

    sns.scatterplot(data=df[df["nbeads"] == 30], x="del_t", y="time", ax=ax, label="30")
    sns.scatterplot(data=df[df["nbeads"] == 40], x="del_t", y="time", ax=ax, label="40")

    fig.savefig(f'{name}_time_test.png')


if __name__ == '__main__':
    time_test('time_test')




