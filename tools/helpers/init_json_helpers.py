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


def run_layer(common_data, in_queue, out_queue, seed_increment, WL_sbias, exe):
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
        common_data['config_in'], common_data['config_out'] = int(format_bits(bits_in), 2), int(format_bits(bits_out),
                                                                                                2)

        # Set the seed
        common_data['seeds'] = [count, seed_increment]

        # Get sbias value for native index from wang-landau (WL) trajectory run along with rc values at different
        # quartiles for this run
        sbias, rc1, rc2, rc3 = wang_landau(common_data, input_hdf5, output_name)

        # Optional staircase potential, initially the bond pair (bp) is not staircased
        stair_bp = None
        stair_rc_list = None

        # if sbias from WL test larger than a threshold value WL_sbias then use staircase potential for this bond
        if sbias > WL_sbias:
            stair_bp = common_data['nonlocal_bonds'][i]
            print(f"Do Staircase on {stair_bp}")
            # Process rc values for staircase from prior WL run; ensure values are in descending order
            stair_rc_list = sorted([round(el, 2) for el in (rc3, rc2, rc1)], reverse=1)

        # optional staircase potential
        if stair_bp:
            data = copy.deepcopy(common_data)
            # iterate through the different staircased rc values
            for j in range(len(stair_rc_list)):

                if j > 0:
                    data['stair'] = stair_rc_list[j - 1]

                if stair_rc_list[j] != stair_rc_list[-1]:
                    output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}_{stair_rc_list[j]}'
                else:
                    output_name = f'hybridmc_{layer}_{format_bits(bits_in)}_{format_bits(bits_out)}'

                data['rc'] = stair_rc_list[j]
                run_sim(data, input_hdf5, output_name, exe)
                input_hdf5 = f'{output_name}.h5'

        else:
            run_sim(common_data, input_hdf5, output_name, exe)

        out_queue.put((bits_out, f'{output_name}.h5'))
