#!/usr/bin/env python3
# runTransition.py runs specified transition to compute dS and the mfpts
#
# example of how to run:
# python ../tools/runTransition.py --json input.json
#
# note that run.py creates a new dir which it enters to generate the output

import argparse
from hybridmc.py_tools.helpers.data_processing_helpers import *
from hybridmc.py_tools.helpers.run_layer_helpers import run_sim, run_stairs
from hybridmc.py_tools.post_processing import get_error, diff_s_bias, mfpt
import multiprocessing as mp


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
    os.chdir('../../')


def get_transition_data(output_name, input_name, exe_name, data, input_hdf5, stair_bp, stair_rc_list, i):
    iter_output_name = f'{output_name}_{i}'
    if os.path.isdir(iter_output_name):
        print(f'{iter_output_name} already exists; saved as old version with given version code')
        os.rename(src=iter_output_name, dst=f"{iter_output_name}.old")

    # make new directory for output for this iteration for the transition simulation
    os.mkdir(iter_output_name)
    os.chdir(iter_output_name)

    # since we are moving into new directory within the loop, change input hdf5 path and exe accordingly
    if input_hdf5 is not None:
        input_hdf5 = f'../{input_name}'

    exe_name = f'../{exe_name}'

    data['seeds'] = [seed + i for seed in data['seeds']]

    run_TransitionProcess(exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list)


# create class to call for multiprocessing purposes
class OutputDataMulti:
    def __init__(self, exe_name, data, input_name, output_name, stair_bp, stair_rc_list):
        self.exe_name = exe_name
        self.data = data
        self.input_hdf5 = input_name
        self.output_name = output_name
        self.stair_bp = stair_bp
        self.stair_rc_list = stair_rc_list

    def __call__(self, i):
        return get_transition_data(self.output_name, self.input_hdf5, self.exe_name, self.data, self.input_hdf5,
                                   self.stair_bp,
                                   self.stair_rc_list, i)


def run_TransitionMulti(json_name, input_name, exe_name, n_iterations=2):
    input_json_name = os.path.basename(json_name)
    data = extract_json(input_json_name)

    output_name = input_json_name.replace('.json', '')

    data = process_json(data, output_name)

    # input hdf5 is None if zeroth layer transition
    if output_name.split('_')[1] == '0':
        input_hdf5 = None

    # otherwise set the input hdf5
    else:
        # input_hdf5 = f'{output_name}.h5'
        input_hdf5 = input_name

    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)

    output_data_multi = OutputDataMulti(exe_name, data, input_hdf5, output_name, stair_bp, stair_rc_list)
    with mp.Pool() as pool:
        pool.map(output_data_multi, range(n_iterations))

    print('Done Running all the simulation iterations')


def extract_json(input_json_name):
    with open(input_json_name, 'r') as input_json:
        data = json.load(input_json)
    return data


def process_json(input_data, output_name):
    # obtain the bits in and bits out strings for the transition
    bitstring_in = output_name.split('_')[2]
    bitstring_out = output_name.split('_')[3]
    input_data["permanent_bonds"] = bonds_from_bitstring(bitstring_in, input_data["nonlocal_bonds"])

    input_data["transient_bonds"] = bonds_from_bitstring(bitstring_subtract(bitstring_out, bitstring_in),
                                                         input_data["nonlocal_bonds"])

    input_data['config_in'] = int(bitstring_in, 2)
    input_data['config_out'] = int(bitstring_out, 2)

    return input_data


def run_TransitionSerial(json_name, input_name, exe_name, n_iterations=2):
    input_json_name = os.path.basename(json_name)
    data = extract_json(input_json_name)

    output_name = input_json_name.replace('.json', '')  # input hdf5 is None if zeroth layer transition

    if output_name.split('_')[1] == '0':
        input_hdf5 = None

    # otherwise set the input hdf5
    else:
        # input_hdf5 = f'{output_name}.h5'
        input_hdf5 = input_name
        # input_hdf5 = f'{output_name}.h5'

    stair_bp, stair_rc_list = stair_check(data, output_name, input_hdf5)

    # loop here with all the different seeds
    for i in range(n_iterations):
        get_transition_data(output_name, exe_name, data, input_hdf5, stair_bp, stair_rc_list, i)


def main(args):
    # obtain mfpt and sbias value by running MD simulation for the transition
    run_TransitionMulti(json_name=args.json, input_name=args.input, exe_name=args.exe, n_iterations=args.n_iterations)

    # run_TransitionSerial(json_name=args.json, exe_name=args.exe)

    # set some more default values for parameters in case they are not provided in json file
    data = extract_json(args.json)
    set_defaults(data,
                 defaults={
                     "D_mean": 0.0394,
                     "D_std_perc": 1,
                     "eps": 3
                 })

    # compile all data into a summary csv file
    get_error.summarize_transition(D_mean=data["D_mean"], D_std_perc=data["D_std_perc"], eps=data["eps"])


if __name__ == '__main__':
    from shutil import which
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='master h5 input file', default='hybridmc_3_0011100000_0011110000.h5')
    parser.add_argument('--json', help='master json input file', default='hybridmc_4_0011110000_0011110010.json')
    parser.add_argument('--exe', help='hybridmc executable', default=which("hybridmc.exe"))

    parser.add_argument('--n_iterations', help='number of iterations to run the transition for averaging',
                        type=int, default=4)

    args = parser.parse_args()

    main(args)
