import copy
import json
import os
import subprocess
import sys
from numpy import searchsorted
from .mfpt_helpers import get_dists
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

    # Exit if this trajectory has already been run. Just another check for safety.
    if os.path.exists(hdf5_name):
        print(f"Simulation results for transition {output_name} already generated")
        return

    # Create json file input for the HMC program to run this trajectory.
    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # Obtain real path for the HMC executable in case it is not already
    exe = os.path.realpath(exe)

    # for layer = 1 or greater,
    command = [exe, json_name, hdf5_name]
    if input_hdf5 is not None:
        command += ['--input-file', input_hdf5]

    print([os.path.relpath(el) for el in command])
    sys.stdout.flush()

    log_name = '{}.log'.format(output_name)
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                       stderr=subprocess.STDOUT)

    sys.stdout.flush()


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
    t_bonds = data['transient_bonds']
    nl_bonds = data['nonlocal_bonds']
    min_rc_percentile = data["rc_target_min_percentile"]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    # initialize rc at the largest staircase rc value;
    # and set the target rc to be the smallest staircase rc value
    rc, rc_target = stair_rc_list[0], stair_rc_list[-1]

    # while the simulation has not converged
    converged = False
    while not converged:

        # append rc to output name
        stair_output_name = f'{output_name}_{rc}'

        # update the rc value for the bonds in the transient_bonds list
        update_rc(data, rc)

        # run the simulation for this staircase step
        run_sim(data, input_hdf5, stair_output_name, exe)

        # setup input file for next step with new rc: it is the output file from the previous step
        input_hdf5 = f'{stair_output_name}.h5'
        # set the outer rc for the next step to be the current rc
        data["stair"] = rc

        # Convergence Check:

        # get the list of distances for the last simulation
        dist_t_active = sorted(get_dists(stair_output_name, t_ind))

        # Figure out what percentile of this distance vector the rc target lies
        rc_target_percentile = searchsorted(dist_t_active, rc_target) / len(dist_t_active)

        # simulation has converged if the rc target percentile is small enough: outer-wall close enough to rc
        if min_rc_percentile <= rc_target_percentile:
            converged = True
            # run last sim with final target rc
            rc = rc_target

        # otherwise keep pushing outer-wall in
        else:
            # find the rc that corresponds to the min_rc_percentile of the distance: set as the new rc
            rc = round(dist_t_active[int(min_rc_percentile * len(dist_t_active))], 1)

            # exit to final sim if it is the target rc
            if rc <= rc_target:
                converged = True

            # Go to start of loop again

    # do final run with rc = target rc
    # update the rc value for the bonds in the transient_bonds list
    update_rc(data, rc)
    run_sim(data, input_hdf5, output_name, exe)
