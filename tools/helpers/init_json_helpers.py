import copy
import json
import os
import HMC
import subprocess
from .data_processing_helpers import *


def wang_landau(nonlocal_bonds_i, data, seed_increment, input_hdf5,
                output_name, bits_in, bits_out, bonds_in, count):
    data['transient_bonds'] = [nonlocal_bonds_i]
    data['permanent_bonds'] = bonds_in
    data['config_in'] = int(format_bits(bits_in), 2)
    data['config_out'] = int(format_bits(bits_out), 2)
    data['seeds'] = [count, seed_increment]

    json_name = f'{output_name}.json'
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # for layer = 1 or greater,
    command = ["wang_landau", json_name]
    json_name = os.path.realpath(json_name)
    input_name = None

    if input_hdf5:
        command += ['--input-file', input_hdf5]
        input_name = os.path.realpath(input_hdf5)

    print(f"Running Wang-Landau procees for {json_name}")

    return HMC.WL_process(json_name, input_name)


def run_sim(nonlocal_bonds_i, data, seed_increment, input_hdf5,
            output_name, bits_in, bits_out, bonds_in, count, exe):
    data['transient_bonds'] = [nonlocal_bonds_i]
    data['permanent_bonds'] = bonds_in
    data['config_in'] = int(format_bits(bits_in), 2)
    data['config_out'] = int(format_bits(bits_out), 2)
    data['seeds'] = [count, seed_increment]

    hdf5_name = f'{output_name}.h5'

    if os.path.exists(hdf5_name):
        print(f"Simulation results for {data['config_in']} to {data['config_out']} transition already generated")
        return

    json_name = f'{output_name}.json'
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # for layer = 1 or greater,
    exe = os.path.realpath(exe)
    json_name = os.path.realpath(json_name)
    hdf5_name = os.path.realpath(hdf5_name)

    # for layer = 1 or greater,
    command = [exe, json_name, hdf5_name]
    if input_hdf5 is not None:
        command += ['--input-file', input_hdf5]

    print(command)
    sys.stdout.flush()

    log_name = '{}.log'.format(output_name)
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                       stderr=subprocess.STDOUT)


def run_layer(common_data, in_queue, out_queue, seed_increment, WL_sbias):
    while True:
        item = in_queue.get()
        if item is None:
            break

        i, bits_in, count, input_hdf5 = item

        # copy of initial state for every nonbonded pair of beads
        bits_out = list(bits_in)

        # flip bit to form bond
        bits_out[i] = True
        bits_out = tuple(bits_out)
        bonds_in = bits_to_bonds(bits_in, common_data['nonlocal_bonds'])
        layer = len(bonds_in)

        # optional staircase potential, initially not on
        bp = None

        output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}'

        sbias_0, sbias_1, rc1, rc2, rc3 = wang_landau(common_data['nonlocal_bonds'][i], common_data, seed_increment,
                                                      input_hdf5, output_name, bits_in, bits_out,
                                                      bonds_in, count)

        # obtain rc values for staircasing; ensure values are in descending order
        rc = sorted([round(el, 2) for el in (rc3, rc2, rc1)], reverse=1)
        sbias = sbias_0 - sbias_1

        # if sbias from wang landau test larger than a threshold value WL_sbias then staircase potential for this bond
        if sbias > WL_sbias:
            bp = nonlocal_bonds[i]
            print(f"Do Staircase on {bp}")

        # optional staircase potential
        if bp:
            data = copy.deepcopy(common_data)
            # iterate through the different staircased rc values
            for j in range(len(rc)):

                if j > 0:
                    data['stair'] = rc[j - 1]

                if rc[j] is not rc[-1]:
                    output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}_{rc[j]}'
                else:
                    output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}'

                data['rc'] = rc[j]
                run_sim(nonlocal_bonds[i], data, seed_increment,
                        input_hdf5, output_name, bits_in, bits_out,
                        bonds_in, count, sbias_1)
                input_hdf5 = f'{output_name}.h5'

        else:
            run_sim(nonlocal_bonds[i], common_data, seed_increment,
                    input_hdf5, output_name, bits_in, bits_out,
                    bonds_in, count, sbias_1)

        out_queue.put((bits_out, f'{output_name}.h5'))
