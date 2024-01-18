#!/usr/bin/env python3
# runTransition.py runs specified transition to compute dS and the mfpts
#
# example of how to run:
# python ../tools/runTransition.py --json input.json
#
# note that run.py creates a new dir which it enters to generate the output
# files, thus, the path of the json input file, init_json.py and the executable
# must be from this newly created dir, and not from the dir from which the
# simulation is initiated


import argparse
import os
import json
from helpers.data_processing_helpers import *
from helpers.run_layer_helpers import run_sim, run_stairs
from post_processing import diff_s_bias, avg_s_bias, mfpt


def post_processing():
    # Obtain the differences in the sbias for each transition
    diff_s_bias.get_diff_sbias()

    # Obtain the average sbias for each bonding state
    avg_s_bias.get_avg_sbias(diff_sbias_csv="diff_s_bias.csv", structure_sim_json=args.json)

    # Obtain the mfpt for each bonding state
    mfpt.get_mfpt()

    # put together the mfpts in one file
    mfpt.compile_mfpts()

def main(args):
    n_iterations = 2
    input_json_name = os.path.basename(args.json)
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)

    print(data)
    output_name = input_json_name.strip('.json')
    print('output name is ', output_name)
    input_hdf5 =  f'{output_name}.h5'
    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)
    print('stair_rc_list is ', stair_rc_list)

    if stair_bp and stair_rc_list:
        run_stairs(data, input_hdf5, args.exe, stair_rc_list)
    else:
        for i in range(n_iterations):
            output_name_i = f'{output_name}.{i}'
            run_sim(data, input_hdf5, output_name_i, args.exe)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='master json input file', default='test.json')
    parser.add_argument('--exe', help='hybridmc executable', default="../../../release/hybridmc")
    args = parser.parse_args()

    #main(args)
    post_processing()
