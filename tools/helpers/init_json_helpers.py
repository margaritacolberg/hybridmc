import copy
import subprocess
from .data_processing_helpers import *


def run_sim(data, input_hdf5, output_name, exe):
    # Set the name for the hdf5 and json files generated
    hdf5_name, json_name = os.path.realpath(f'{output_name}.h5'), os.path.realpath(f'{output_name}.json')

    # Exit if this trajectory has already been run
    if os.path.exists(hdf5_name):
        print(f"Simulation results for {data['config_in']} to {data['config_out']} transition already generated")
        return

    # Create json file input for the HMC program to run this trajectory
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # Obtain real path for the HMC executable
    exe = os.path.realpath(exe)

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


def run_layer(common_data, in_queue, out_queue, seed_increment, exe):
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

        # Define the output files identifier for this transition (used to name hdf5 files, json files and log files)
        output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}'

        # Set the transient bond to form this trajectory
        common_data['transient_bonds'] = [common_data['nonlocal_bonds'][i]]

        # Set the existing bonds as permanent bonds
        common_data['permanent_bonds'] = bonds_in

        # Obtain configuration information
        common_data['config_in'] = int(format_bits(bits_in), 2)
        common_data['config_out'] = int(format_bits(bits_out), 2)

        # Set the seed
        common_data['seeds'] = [count, seed_increment]

        # Check if staircase needed and do a staircase run for structure; if not do regular run
        if not run_stairs(common_data, exe, input_hdf5, output_name):
            run_sim(common_data, input_hdf5, output_name, exe)

        out_queue.put((bits_out, f'{output_name}.h5'))


def run_stairs(common_data, exe, input_hdf5, output_name):

    # Check for stairs using wang_landau (WL) run
    stair_bp, stair_rc_list = stair_check(common_data, output_name, input_hdf5)

    # optional staircase potential
    if stair_bp:
        data = copy.deepcopy(common_data)
        # iterate through the different staircased rc values
        for j in range(len(stair_rc_list)):

            if j > 0:
                data['stair'] = stair_rc_list[j - 1]

            if stair_rc_list[j] != stair_rc_list[-1]:
                stair_output_name = f'{output_name}_{stair_rc_list[j]}'
            else:
                stair_output_name = output_name

            data['rc'] = stair_rc_list[j]
            run_sim(data, input_hdf5, stair_output_name, exe)
            input_hdf5 = f'{stair_output_name}.h5'

        print(f"Finished staircase run for {common_data['transient_bonds']}")
        return True

    else:
        print(f"No staircase run for {common_data['transient_bonds']}")
        return False
