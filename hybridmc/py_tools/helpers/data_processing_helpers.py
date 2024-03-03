import sys
import os
import json

import h5py
import wang_landau as WL


def n_sims(nbonds):
    return nbonds * 2 ** (nbonds - 1)


def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data


def make_rc_tuple(nonlocal_bonds, rc):
    """
    Function to produce nonlocal bonds list with each nonlocal bond list element containing the rc as the third
    element

    Parameters
    ----------
    nonlocal_bonds -- List of Lists containing the indices of beads bonded to each other, for example [ [1, 9], [4, 15]].
    These elements may potentially have the bond length (rc) as the third element as well or not,
    for example [[1, 9, 2.3], [2, 6]].

    rc -- Default rc to be appended to each nonlocal_bonds element in case they do not have this information

    Returns
    -------
    nonlocal_bonds -- The same list but with each list element having the bond length (rc) as their third element.
    """

    # loop through each bonded bead pair
    for el in nonlocal_bonds:

        # check if rc for the bond is already given
        if len(el) != 3:
            # if not then append rc
            el.append(rc)

    return nonlocal_bonds


def sort_triplet(bond_list):
    """
swap bond i and j in (i, j, rc) triplets within bond list if i > j

Parameters
----------
bond_list -- List of Lists containing the indices of beads bonded to each other, for example [ [1, 9], [4, 15]]. These
elements may potentially have the bond length (rc) as the third element as well or not, for example [[1, 9, 2.3], [2, 6]].

"""

    for el in bond_list:
        if el[0] > el[1]:
            el[0], el[1] = el[1], el[0]

    return bond_list


def stairs_rcval_check(bond_pair, stair_rcs, nonlocal_bonds):
    """
    function to check if staircase distances are appropriate

    """
    rc_check = 0
    for el in nonlocal_bonds:
        if el[:-1] == bond_pair:
            rc_check = el[-1]
    problem_rc = []
    for el in stair_rcs:
        if el < rc_check:
            problem_rc.append(el)
    if problem_rc:
        sys.exit((
            f"Value error: rc values given ({problem_rc}) are too small, no rc can be less than smallest allowed ({rc_check})"))


def bits_to_bonds(bits, nonlocal_bonds):
    """
    function to convert bits to bonds
    Parameters
    ----------
    bits: iterable of bool
    nonlocal_bonds: List[List]: List of all the possible nonlocal bonds

    Returns
    -------
    List of bits pattern described bonds

    """
    # obtain index of the turned on bonds (postion of every bits = 1)
    bits = [i for i, bit in enumerate(bits) if bit]

    # Filter out bonds with indices in bits from nonlocal bonds (these are on)
    bonds = [nonlocal_bonds[i] for i in bits]
    return bonds


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))


def bonds_from_bitstring(bitstring, nonlocal_bonds):
    """
    Uses bitstring as a mask to extract relevant bonds from the list of all bonds in nonlocal bonds

    Example: if bitstring = '001', nonlocal_bonds = [[1,2], [3,4], [5,6]]
    result = [[5,6]]

    Parameters
    ----------
    bitstring: str: String specifying mask to extract relevant bonds from in nonlocal bonds
    nonlocal_bonds: list: List of all non-local bonds

    Returns
    -------
    Extracted all bonds specified by bitstring in a list

    """

    return bits_to_bonds((int(el) for el in bitstring), nonlocal_bonds)


def bitstring_subtract(bitstring_out, bitstring_in):
    """
    Subtracts the bitstrings given (out - in). It is an element by element subtraction.
    Each string in the bistrings represents a bit (0 or 1)

    Parameters
    ----------
    bitstring_out
    bitstring_in

    Returns
    -------
    The subtracted bitstring
    """
    result = ''
    for i in range(len(bitstring_out)):
        result += str(int(bitstring_out[i]) - int(bitstring_in[i]))

    return result


def extract_sbias(file_path):
    with h5py.File(file_path, 'r') as f:
        return f['s_bias'][0]


def set_defaults(data, defaults):
    """
    Function to set default key: values pairs to the data dictionary in case the key does not exist.
    Parameters
    ----------
    data: Dict: The main data directory
    defaults: Dict: Contains default key: value pairs to add to data.

    Returns
    -------
    None
    """
    for el in defaults.items():
        data.setdefault(*el)


def stair_check(data, output_name, input_hdf5):
    """
    Function to check if staircase potential is needed for a given bond pair.
    This is done by running a Wang-Landau process on the bond pair and checking if
    the sbias value is larger than a threshold value. If so, then staircase potential is needed.

    Parameters
    ----------
    data: dict -- dictionary containing the input parameters for the wang-landau run
    output_name: str -- name of the output json file
    input_hdf5: str -- name of the input hdf5 file

    Returns
    -------
    stair_bp: list -- list containing the bond pair to use a staircase potential on
    stair_rc_list: list -- list containing the rc values for the staircase potential
    """

    hdf5_name, json_name = os.path.realpath(f'{output_name}.h5'), os.path.realpath(f'{output_name}.json')
    print(f"Running Wang-Landau process for {os.path.relpath(json_name)} to check if staircase potential needed")

    # find what configurational change is happening for this simulation
    config_in, config_out = json_name.replace('.json', '').split('_')[-2], json_name.strip('.json').split('_')[-1]

    # Create json input for the wang_landau run

    with open(json_name, 'w') as output_json:
        json.dump(data, output_json)

    # Get sbias value for native index from wang-landau (WL) trajectory run along with rc values at different
    # quartiles for this run

    # find pre-existing stair sims rc values either in h5 or tmp files
    existing_stairs = (if_stair(output_name, os.listdir(), end_tag='.h5') +
                       if_stair(output_name, os.listdir(), end_tag='.h5.tmp'))

    if len(existing_stairs) > 0:
        sbias = data["WL_sbias"]
        rc1, rc2 = sorted(existing_stairs)[0], data["transient_bonds"][-1][-1]
        rc3 = 0.5 * (rc1 + rc2)

    else:
        sbias, rc1, rc2, rc3 = WL.WL_process(json_name, input_hdf5)

    # Optional staircase potential, initially the bond pair (bp) is not staircased
    stair_bp = None
    stair_rc_list = None

    # if sbias from WL test larger than or equal to threshold value WL_sbias then use staircase potential for this bond
    if sbias >= data['WL_sbias']:
        stair_bp = data['transient_bonds']
        print(f"Do Staircase on {stair_bp} because sbias is {sbias} for transition {config_in} to {config_out}")
        # Process rc values for staircase from prior WL run; ensure values are in descending order
        stair_rc_list = sorted([round(el, 2) for el in (rc3, rc2, rc1)], reverse=1)
    else:
        print(f"NO Staircase because sbias is {sbias} for transition {config_in} to {config_out}")

    return stair_bp, stair_rc_list


def update_rc(data, new_rc):
    """
    Update the rc value for the bonds in the transient_bonds list and nonlocal_bonds list

    Parameters
    ----------
    data: dict: JSON input for the HMC program as a python dictionary
    new_rc: float: new rc value to be set for the bonds in the transient_bonds list

    Returns
    -------
    None
    """
    # update the rc value for the bonds in the transient_bonds list
    for bonds in data['nonlocal_bonds']:
        # check if this bond is in the transient_bonds list
        if bonds in data['transient_bonds']:
            # if yes then update the rc value
            bonds[-1] = new_rc
            data['transient_bonds'] = [bonds]


def get_mcMoves(data: dict):
    """
    Function to determine how many monte-carlo moves to make based on the number of bonding restrictions present

    Parameters
    ----------
    data: dict: JSON input for the HMC program as a python dictionary

    Returns
    -------
    mcMoves: int: number of monte-carlo moves to make
    """
    return round(1000 * 1 - (len(data['permanent_bonds']) / len(data['nonlocal_bonds'])))


def if_stair(ref_sim_id, files, end_tag='.json'):
    """
    Function to check if the file_path has associated staircased steps simulations results. If yes, then the paths of these
    intermediate steps' csv files with their mfpt information are compiled into a list and returned.

    :param file_path: the final step mfpt simulation id
    :param files: the list of files in directory of interest
    :return: list[float] containing all the intermediate staircase rc values if they exist
    """

    # initialize output list merging all intermediate stair mfpt csv paths
    stair_rc_list = []

    # loop through json files in directory
    for file in files:
        if file.endswith(end_tag):
            # obtain simulation tag for this file
            sim_id = file.replace(end_tag, '')
            # if the reference tag and this are the same then add the filepath to the output list
            if ref_sim_id.split('_') == sim_id.split('_')[:-1]:
                stair_rc_list.append(float(sim_id.split("_")[-1]))

    return stair_rc_list
