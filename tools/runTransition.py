#!/usr/bin/env python3
# runTransition.py runs specified transition to compute dS and the mfpts
#
# example of how to run:
# python ../tools/runTransition.py --json input.json
#
# note that run.py creates a new dir which it enters to generate the output

import argparse
import os
import json
from helpers.data_processing_helpers import *
from helpers.run_layer_helpers import run_sim, run_stairs
from post_processing import diff_s_bias, avg_s_bias, mfpt
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import glob
from pandas import read_csv
import csv
import numpy as np


def run_Transition(args):
    n_iterations = 2
    input_json_name = os.path.basename(args.json)
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)

    output_name = input_json_name.strip('.json')
    input_hdf5 = f'{output_name}.h5'
    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)
    #stair_bp, stair_rc_list = 0, [0]

    # since we are moving into new directory within the loop, change input hdf5 path and exe accordingly
    input_hdf5 = f'../{input_hdf5}'
    args.exe = f'../{args.exe}'

    # loop here with all the different seeds
    for i in range(n_iterations):
        iter_output_name = f'{output_name}_{i}'
        if os.path.isdir(iter_output_name):
            print(f'{iter_output_name} already exists; saved as old version with given version code')
            os.rename(src=iter_output_name, dst=f"{iter_output_name}.old")

        # make new directory for output for this iteration for the transition simulation
        os.mkdir(iter_output_name)
        os.chdir(iter_output_name)

        data['seeds'] = [i + 1 for i in data['seeds']]


        #exe.map(run_TransitionProcess, (args, data, input_hdf5, output_name, stair_bp, stair_rc_list))

        run_TransitionProcess(args, data, input_hdf5, output_name, stair_bp, stair_rc_list)


def run_TransitionProcess(args, data, input_hdf5, output_name, stair_bp, stair_rc_list):
    if stair_bp and stair_rc_list:
        run_stairs(data, input_hdf5, output_name, args.exe, stair_rc_list)

    else:
        run_sim(data, input_hdf5, output_name, args.exe)

    post_processing()
    os.chdir('../')


def post_processing():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias()

    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt()

    # put together the mfpts in one file
    mfpt.compile_mfpts()

def compute_averages():
    src_bias, src_mfpt = '*/diff_s_bias.csv', '*/mfpt.csv'

    diff_sbias = []
    inner_mfpts = []
    outer_mfpts = []
    # average diff_s_bias
    for csv_file in glob.glob(src_bias):
        sbias_data = read_csv(csv_file, header=None)
        diff_sbias.append(sbias_data.iloc[0, 2])

    for csv_file in glob.glob(src_mfpt):
        mfpts_data = read_csv(csv_file)
        inner_mfpts.append(mfpts_data.iloc[:, 3].values[0])
        outer_mfpts.append(mfpts_data.iloc[:, 4].values[0])

    output = []
    for data in [diff_sbias, inner_mfpts, outer_mfpts]:
        mean, var = np.mean(data), np.var(data)
        percent_rel_e = 100 * var / mean
        output.append([mean, var, percent_rel_e])


    csv_name = 'summary_data.csv'
    with open(csv_name, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='hybridmc_3_1000000110_1010000110.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../../../release/hybridmc")
    args = parser.parse_args()

    run_Transition(args)
