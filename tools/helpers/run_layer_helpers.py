import copy
import json
import os
import csv
import subprocess
import sys
import h5py
import numpy as np
from .data_processing_helpers import update_rc
from .mfpt_helpers import find_knots 
from .mfpt_helpers import fpt
from .mfpt_helpers import fpt_per_bead_pair


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

    print([os.path.relpath(el) for el in command])
    sys.stdout.flush()

    log_name = '{}.log'.format(output_name)
    with open(log_name, 'w') as output_log:
        subprocess.run(command, check=True, stdout=output_log,
                       stderr=subprocess.STDOUT)

    sys.stdout.flush()

def inner_stair_fpt(json_in,hdf5_in):
    with open(json_in, 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']

    t_bonds = data['transient_bonds']
    p_bonds = data['permanent_bonds']
    nl_bonds = data['nonlocal_bonds']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    if p_bonds:
        p_ind = [i for i, bond in enumerate(nl_bonds) for j in
                 range(len(p_bonds)) if bond == p_bonds[j]]
    else:
        p_ind = []

    with h5py.File(hdf5_in, 'r') as f:
        dist = f['dist'][:]

    dist = dist.T
    nbonds = len(nl_bonds)

    t_on = []

    for j in range(len(dist[t_ind])):
        if dist[t_ind][j] < rc_transient:
            t_on.append(dist[t_ind][j])

    nknots = 12
    t_on = np.array(t_on)
    fpt_on = fpt_per_bead_pair(t_on, nknots, rh, rc_transient, True)[0]
    return t_ind, fpt_on


def knots_wrap(json_in, hdf5_in, rc_max):
    with open(json_in, 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']

    t_bonds = data['transient_bonds']
    p_bonds = data['permanent_bonds']
    nl_bonds = data['nonlocal_bonds']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    if p_bonds:
        p_ind = [i for i, bond in enumerate(nl_bonds) for j in
                 range(len(p_bonds)) if bond == p_bonds[j]]
    else:
        p_ind = []

    with h5py.File(hdf5_in, 'r') as f:
        dist = f['dist'][:]


    dist = dist.T
    nbonds = len(nl_bonds)

    t_on = []
    t_off = []

    for j in range(len(dist[t_ind])):
        if dist[t_ind][j] < rc_transient:
            t_on.append(dist[t_ind][j])
        else:
            t_off.append(dist[t_ind][j])

    t_off = np.array(t_off)

    if rc_max == -1.:
        rc_max = np.max(t_off)

    return find_knots(t_off, rc_transient, rc_max)


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
        input_hdf5 = f'{stair_output_name}.h5'


    xknots_total = []
    yknots_total = []
    for j in range(len(stair_rc_list)):

        if j > 0:
            data['stair'] = stair_rc_list[j - 1]

        # So long as the rc is not at the target we append the intermediate rc to the output name
        if stair_rc_list[j] != stair_rc_list[-1]:
            stair_output_name = f'{output_name}_{stair_rc_list[j]}'
        else:
            stair_output_name = output_name

        input_hdf5 = f'{stair_output_name}.h5'
        json_name = f'{stair_output_name}.json'

        if j == 0:
            rc_max = -1.
        else:
            rc_max = stair_rc_list[j-1]

        nknots, x_knot, y_knot = knots_wrap(json_name, input_hdf5,rc_max) 

        if j == 0:
            xknots_total = x_knot.copy()
            yknots_total = y_knot.copy()
        else:
            # exclude last x_knot, y_knot
            l = len(x_knot)
            last_y = yknots_total[0]
            for i in range(l-2,-1,-1):
                xknots_total = np.insert(xknots_total, 0, x_knot[i])
                yknots_total = np.insert(yknots_total, 0, y_knot[i] + last_y)

    t_ind, inner_fpt = inner_stair_fpt(json_name, input_hdf5) 
    outer_fpt = fpt(xknots_total, yknots_total, xknots_total[0], xknots_total[-1], False)[0]

    output = []
    output.append([t_ind, inner_fpt,outer_fpt])
    csv_name = f'{output_name}.csv'
    with open(csv_name, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)

    return
