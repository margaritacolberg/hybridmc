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


def run_Transition(args):
    n_iterations = 2
    input_json_name = os.path.basename(args.json)
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)

    output_name = input_json_name.strip('.json')
    input_hdf5 = f'{output_name}.h5'
    #stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)
    stair_bp, stair_rc_list = 0, [0]

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
        run_stairs(data, input_hdf5, args.exe, stair_rc_list)

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='hybridmc_3_0000011100_0010011100.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../../../release/hybridmc")
    args = parser.parse_args()

    run_Transition(args)
