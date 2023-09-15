import copy
import json
import os
import subprocess
import sys
from .data_processing_helpers import update_rc


def run_sim(data, input_hdf5, output_name, exe):
    """
    Run a single simulation using the HMC program for one transition in one layer
    Parameters
    ----------
    data: dict: JSON input for the HMC program as a python dictionary
    input_hdf5: str: name of the input hdf5 file
    output_name: str: name of the output file
    exe: str: path to the HMC executable

    Returns
    -------
    None
    """

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


def run_stairs(common_data, input_hdf5, output_name, exe, stair_rc_list):
    """
    Run staircase simulations for a given transition:

    Use rc values in stair_rc_list to define staircase potential for this transition.
    Pipe each of the steps to run_sim and use results to eventually push beads together within a
    distance of rc -- their target bonding rc.

    This function returns nothing, instead it will produce files for each staircase step. The final step
    output files will look like the normal non-staircase files but the intermediate outputs will have
    the staircase inner rc value appended to the end of the file name.

    Parameters
    ----------
    common_data: dict: JSON input for the HMC program as a python dictionary
    input_hdf5: str: name of the input hdf5 file
    output_name: str: name of the output file
    exe: str: path to the HMC executable
    stair_rc_list: list: list of rc values to define staircase potential for this transition

    Returns
    -------
    None
    """

    data = copy.deepcopy(common_data)
    # iterate through the different staircase rc values
    for j in range(len(stair_rc_list)):

        # Once we push beads within the largest staircase rc, we can set that rc to be outer wall
        # We call this boundary the "stair"
        if j > 0:
            data['stair'] = stair_rc_list[j - 1]

        # So long as the rc is not at the target we append the intermediate rc to the output name
        if stair_rc_list[j] != stair_rc_list[-1]:
            stair_output_name = f'{output_name}_{stair_rc_list[j]}'
        else:
            stair_output_name = output_name

        # update the rc value for the bonds in the transient_bonds list
        update_rc(data, stair_rc_list[j])

        # run the simulation for this staircase step
        run_sim(data, input_hdf5, stair_output_name, exe)

        # update the input hdf5 file for the next step
        input_hdf5 = f'{stair_output_name}.h5'

    print(f"Finished staircase run for {common_data['transient_bonds']}")

    return
