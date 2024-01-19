#!/usr/bin/env python3
# runTransition.py runs specified transition to compute dS and the mfpts
#
# example of how to run:
# python ../tools/runTransition.py --json input.json
#
# note that run.py creates a new dir which it enters to generate the output

import argparse
from helpers.data_processing_helpers import *
from helpers.run_layer_helpers import run_sim, run_stairs
from post_processing import diff_s_bias, mfpt
import multiprocessing as mp
import glob
from pandas import read_csv
from csv import DictWriter
import numpy as np


def post_processing():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias()

    # Obtain the mfpt for each bonding state (only one transition so do serially)
    mfpt.get_mfpt_serial()

    # put together the mfpts in one file
    mfpt.compile_mfpts()


def run_TransitionProcess(exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list):
    if stair_bp and stair_rc_list:
        run_stairs(data, input_hdf5, output_name, exe_name, stair_rc_list)

    else:
        run_sim(data, input_hdf5, output_name, exe_name)

    post_processing()
    os.chdir('../')


def get_transition_data(output_name, exe_name, data, input_hdf5, stair_bp, stair_rc_list, i):
    iter_output_name = f'{output_name}_{i}'
    if os.path.isdir(iter_output_name):
        print(f'{iter_output_name} already exists; saved as old version with given version code')
        os.rename(src=iter_output_name, dst=f"{iter_output_name}.old")

    # make new directory for output for this iteration for the transition simulation
    os.mkdir(iter_output_name)
    os.chdir(iter_output_name)

    data['seeds'] = [i + 1 for i in data['seeds']]

    run_TransitionProcess(exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list)


# create class to call for multiprocessing purposes
class OutputDataMulti(object):
    def __init__(self, exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list):
        self.exe_name = exe_name
        self.data = data
        self.input_hdf5 = input_hdf5
        self.output_name = output_name
        self.stair_bp = stair_bp
        self.stair_rc_list = stair_rc_list

    def __call__(self, i):
        return get_transition_data(self.output_name, self.exe_name, self.data, self.input_hdf5, self.stair_bp,
                                   self.stair_rc_list, i)


def run_TransitionMulti(json_name, exe_name, n_iterations=4):
    input_json_name = os.path.basename(json_name)
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)

    output_name = input_json_name.strip('.json')
    input_hdf5 = f'{output_name}.h5'
    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)
    # stair_bp, stair_rc_list = 0, [0]

    # since we are moving into new directory within the loop, change input hdf5 path and exe accordingly
    input_hdf5 = f'../{input_hdf5}'
    exe_name = f'../{exe_name}'

    output_data_multi = OutputDataMulti(exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list)
    with mp.Pool() as pool:
        pool.map(output_data_multi, range(n_iterations))

    print('Done Running all the simulation iterations')


def run_TransitionSerial(json_name, exe_name, n_iterations=2):
    input_json_name = os.path.basename(json_name)
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)

    output_name = input_json_name.strip('.json')
    input_hdf5 = f'{output_name}.h5'
    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)
    # stair_bp, stair_rc_list = 0, [0]

    # since we are moving into new directory within the loop, change input hdf5 path and exe accordingly
    input_hdf5 = f'../{input_hdf5}'
    exe_name = f'../{exe_name}'

    # loop here with all the different seeds
    for i in range(n_iterations):
        get_transition_data(output_name, exe_name, data, input_hdf5, stair_bp, stair_rc_list, i)


def summarize_transition():
    # define patterns for diff s bias and mfpt files
    src_bias, src_mfpt = '*/diff_s_bias.csv', '*/mfpt.csv'

    # intitialize lists to hold sbias and mfpts for the transition
    diff_sbias = []
    inner_mfpts = []
    outer_mfpts = []

    # average diff_s_bias compilation
    for csv_file in glob.glob(src_bias):
        sbias_data = read_csv(csv_file, header=None)
        diff_sbias.append(sbias_data.iloc[0, 2])

    # mfpt compilation
    for csv_file in glob.glob(src_mfpt):
        mfpts_data = read_csv(csv_file)
        inner_mfpts.append(mfpts_data.iloc[:, 3].values[0])
        outer_mfpts.append(mfpts_data.iloc[:, 4].values[0])

    # put values together in a list of dictionaries
    output = [
        {'ID': 'sbias', 'vals': diff_sbias},
        {'ID': 'inner_mfpts', 'vals': inner_mfpts},
        {'ID': 'outer_mfpts', 'vals': outer_mfpts}
    ]

    # for each dict in output find relevant stats
    for data in output:
        data['count'], data['mean'], data['var'] = len(data['vals']), np.mean(data['vals']), np.var(data['vals'])
        data['percent_rel_e'] = 100 * data['var'] / data['mean']
        data.pop('vals')

    # write the dict with stats to a csv file
    csv_name = 'summary_data.csv'
    with open(csv_name, 'w') as output_csv:
        writer = DictWriter(output_csv, fieldnames=output[0].keys())
        writer.writeheader()
        writer.writerows(output)

    print('Completed compiling and writing transition stats for sbias and mfpt')


def main(args):
    # obtain mfpt and sbias value by running MD simulation for the transition
    run_TransitionMulti(json_name=args.json, exe_name=args.exe)

    # run_TransitionSerial(json_name=args.json, exe_name=args.exe)

    # compile all data into a summary csv file
    summarize_transition()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='hybridmc_3_1000000110_1010000110.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../../../release/hybridmc")
    args = parser.parse_args()

    main(args)
