import copy
import json
import os
import HMC
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

    print(command)
    #sys.stdout.flush()

    return HMC.WL_process(json_name, input_name)


def run_sim(nonlocal_bonds_i, data, seed_increment, input_hdf5, hdf5_name,
            output_name, bits_in, bits_out, bonds_in, count, init_sbias):
    data['transient_bonds'] = [nonlocal_bonds_i]
    data['permanent_bonds'] = bonds_in
    data['config_in'] = int(format_bits(bits_in), 2)
    data['config_out'] = int(format_bits(bits_out), 2)
    data['seeds'] = [count, seed_increment]

    if os.path.exists(hdf5_name):
        return

    json_name = f'{output_name}.json'
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # for layer = 1 or greater,
    command = ["run_simulation", json_name, hdf5_name]
    json_name = os.path.realpath(json_name)
    hdf5_name = os.path.realpath(hdf5_name)

    input_name = None
    if input_hdf5:
        command += ['--input-file', input_hdf5]
        input_name = os.path.realpath(input_hdf5)

    print(command)
    #sys.stdout.flush()

    HMC.adaptive_convergence(json_name, hdf5_name, init_sbias, input_name)


def run_layer(nonlocal_bonds, common_data, in_queue, out_queue, seed_increment, WL_sbias):
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
        bonds_in = bits_to_bonds(bits_in, nonlocal_bonds)
        layer = len(bonds_in)

        # optional staircase potential, initially not on
        bp = None
        # seperate bond pair from rc in the nonlocal_bonds list
        rc_i = nonlocal_bonds[i][-1]  # rc is always the last element
        nonlocal_bonds_i = nonlocal_bonds[i][:-1]  # so exclude last element to get the nonlocal bonds pair

        output_name = 'hybridmc_{}_{}_{}'.format(layer,
                                                 format_bits(bits_in), format_bits(bits_out))

        sbias_0, sbias_1, rc1, rc2, rc3 = wang_landau(nonlocal_bonds[i], common_data, seed_increment,
                                                      input_hdf5, output_name, bits_in, bits_out,
                                                      bonds_in, count)

        rc = [round(el, 2) for el in (rc3, rc2, rc1)]
        sbias = sbias_0 - sbias_1

        # if sbias from wang landau test larger than a threshold value WL_sbias then staircase potential for this bond
        if sbias > WL_sbias:
            bp = nonlocal_bonds[i]

        # optional staircase potential
        if bp:
            for j in range(len(rc)):
                data = copy.deepcopy(common_data)

                # if there is at least one permanent bond,
                if layer > 0:
                    data['p_rc'] = rc_i

                if j > 0:
                    data['stair_bonds'] = [nonlocal_bonds[i]]
                    data['stair'] = rc[j - 1]

                if j < len(rc) - 1:
                    output_name = 'hybridmc_{}_{}_{}_{}'.format(layer,
                                                                format_bits(bits_in), format_bits(bits_out), rc[j])
                    hdf5_name = '{}.h5'.format(output_name)
                else:
                    output_name = 'hybridmc_{}_{}_{}'.format(layer,
                                                             format_bits(bits_in), format_bits(bits_out))
                    hdf5_name = '{}.h5'.format(output_name)

                data['rc'] = rc[j]
                run_sim(nonlocal_bonds[i], data, seed_increment,
                        input_hdf5, hdf5_name, output_name, bits_in, bits_out,
                        bonds_in, count, sbias_1)
                input_hdf5 = hdf5_name

        else:
            output_name = 'hybridmc_{}_{}_{}'.format(layer,
                                                     format_bits(bits_in), format_bits(bits_out))
            hdf5_name = '{}.h5'.format(output_name)
            run_sim(nonlocal_bonds[i], common_data, seed_increment,
                    input_hdf5, hdf5_name, output_name, bits_in, bits_out,
                    bonds_in, count, sbias_1)

        out_queue.put((bits_out, hdf5_name))
