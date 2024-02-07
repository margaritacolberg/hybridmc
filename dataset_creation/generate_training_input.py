# -*- coding: utf-8 -*-
"""
Created on Fri May 12 17:03:10 2021

@author: Vigneshwar Rajesh
"""
import os
from random import randint as ri, sample, choice as ch
import csv
import json
import sys

cols = ["nbeads", "nbonds", "bead_index"]

# A dict of the form dict(nbonds: dict(row type: number of rows  for this type))
num_rows = {
    3: {
        'total': 1875,
        'mixed': 1524,
        'alpha': 39,
        'beta': 312,
    },
    4: {
        'total': 940,
        'mixed': 739,
        'alpha': 45,
        'beta': 156,
    },
    5: {
        'total': 320,
        'mixed': 216,
        'alpha': 51,
        'beta': 53,
    },
    6: {
        'total': 80,
        'mixed': 42,
        'alpha': 19,
        'beta': 19,
    },
    7: {
        'total': 20,
        'mixed': 12,
        'alpha': 4,
        'beta': 4,
    },
    8: {
        'total': 8,
        'mixed': 6,
        'alpha': 1,
        'beta': 1,
    }

}

num_rows_v2 = {
    3: {
        'total': 1875,
        'mixed': 1524,
        'alpha': 32,
        'beta': 319,
    },
    4: {
        'total': 944,
        'mixed': 739,
        'alpha': 35,
        'beta': 170,
    },
    5: {
        'total': 320,
        'mixed': 216,
        'alpha': 45,
        'beta': 59,
    },
    6: {
        'total': 80,
        'mixed': 42,
        'alpha': 19,
        'beta': 19,
    },
    7: {
        'total': 20,
        'mixed': 12,
        'alpha': 4,
        'beta': 4,
    },
    8: {
        'total': 8,
        'mixed': 6,
        'alpha': 1,
        'beta': 1,
    }

}

common = {
    "m": 6.0,
    "sigma_bb": 1.0,
    "near_min": 1.0,
    "near_max": 1.17,
    "nnear_min": 1.4,
    "nnear_max": 1.67,
    "rh": 1.25,
    "rc": 1.5,
    # "nonlocal_bonds": [[4, 12, 1.6], [6, 10, 1.5]],
    # "permanent_bonds": [],
    "nbeads": 20,
    "tries": 10000,
    "length": 25.0,
    "ncell": 5,
    "nsteps": 1,
    "nsteps_eq": 8,
    "del_t": 1.0,
    "nsteps_wl": 5,
    "del_t_wl": 2.0,
    "gamma": 0.1,
    "gamma_f": 0.0005,
    "write_step": 5,
    "seeds": [3, 1],
    "temp": 1.0,
    "mc_moves": 100,
    "total_iter": 1600,
    "total_iter_eq": 100,
    "pos_scale": 0.5,
    "neg_scale": 0.1,
    "sig_level": 0.1,
    "max_nbonds": 1,
    "max_g_test_count": 10000,
    "flip_req": 1.0,
    "WL_sbias": 5.0,
    "fail_max": 1,
    "req_dists": 50000,
    "rc_target_min_percentile": 0.025

}


def format_bonds(nonlocal_bonds):
    """
    Format the list of non-local bonds such that the smaller index comes first for each pair of indices
    describing a bond and that the bond pair with the smaller first index comes first. So
    [[22,3], [11,3], [9,4]] will become [[3, 11], [3, 22], [4, 9]]
    Parameters
    ----------
    nonlocal_bonds (List[List]): List containing the pairs of values describing different bonds

    Returns
    -------
    Ordered list based on nonlocal bonds

    """
    nonlocal_bonds = [sorted(el) for el in nonlocal_bonds]
    nonlocal_bonds.sort()

    return nonlocal_bonds


# For generating alpha helix values
def get_alpha_helix_row(num_bonds):
    """
    Function for generating alpha helix bonds

    Parameters
    ----------
    num_bonds: int. The number of bonds in this structure
    rows: set. Contains the structure row values (nbeads, nbonds, bead_index)

    Returns
    -------
    tuple
    """

    bonds = []
    # Set the number of beads range based on the number of bonds to ensure tails that are not too long but long enough
    # to have all bonds present
    minbeads, maxbeads = num_bonds * 5, num_bonds * 5 + 5

    # Randomly select the number of beads
    nbeads = ri(minbeads, maxbeads)

    # Randomly pick which index to start the alpha helix pattern between 0 and the last index
    # that allows all bonds to be palced ina n alpha helix pattern (4 * num_bonds - 1)
    ind = ri(0, nbeads - (4 * num_bonds) - 1)

    # create bonds based on the alpha helix pattern until num_bonds pairs of values are generated
    while num_bonds and nbeads - ind >= 4:
        bonds.append([ind, ind + 4])
        ind += 4
        num_bonds -= 1

    # return the tuple containing the alpha helix structure row
    return tuple((nbeads, len(bonds), str(format_bonds(bonds))))


# Some helper functions for making beta sheet values
def get_params(nbeads, mg):
    """
    Parameters
    ----------
    nbeads:int
    mg: int. minimum gap between strands

    Returns
    -------
    strandsize, number of strands, starting position for beads : int, int, int
    """
    start = ri(0, nbeads - 5 - mg)
    new_beads = nbeads - start  # new number of available indices

    max_ns = 1 + ((new_beads - 2) // (mg + 2))  # max strands number found based on new_beads and mingap
    nstrands = ri(2, max_ns)  # randomize number of strands. Min is 2

    # maximum strand size based on nstrands, new_beads. Strandsize Cannot go over 8 beads
    max_sz = min(((new_beads - mg * (nstrands - 1)) // nstrands, 8))
    strandsize = ri(2, max_sz)  # randomize strandsize. Min is 2.

    return strandsize, nstrands, start


def get_strands(nbeads, mg):
    """
    Parameters
    ----------
    nbeads: int
    mg: int. minimum gap between strands

    Returns
    -------
    list of strands indices, number of strands (i.e len(strands)) > 2
    """
    strandsize, nstrands, start = get_params(nbeads, mg)
    strands = []
    for i in range(nstrands):
        # Get the strand indices list and compile them in strands list
        strand = [el for el in range(start, start + strandsize)]
        strands.append(strand)

        # Get start position for the next strand from gaptn and strandsize
        maxgap = nbeads - strand[-1] - 1 - strandsize * (nstrands - i - 1) - mg * (nstrands - i - 2)
        gaptn = ri(mg, maxgap)  # randomize the gap till the next strand.
        start += strandsize + gaptn

    return strands


def get_beta_row(num_bonds, mg=2):
    """
    Function for generating the beta sheet values

    Parameters
    ----------
    num_bonds: int. The number of bonds in this structure
    rows: set. Contains the structure row values (nbeads, nbonds, bead_index)
    mg: int. minimum gap between strands

    Returns
    -------
    tuple
    """
    bonds = []  # initialize bonds list for compiling bond indices
    minbeads, maxbeads = num_bonds * 2 + 2, num_bonds * 5 + 5
    while not bonds:
        nbeads = ri(minbeads, maxbeads)
        strands = get_strands(nbeads, mg)  # for this sheet, assumed to be >= 2
        strand_size = len(strands[1])  # strand size
        sense = ch(['p', 'ap'])  # The sense can be 'p' or 'ap' for parallel of  anti-parallel

        for i in range(len(strands) - 1):
            for j in range(strand_size):
                if len(bonds) != num_bonds:
                    if sense == 'ap':
                        bonds.append([strands[i][j], strands[i + 1][strand_size - 1 - j]])
                    elif sense == 'p':
                        bonds.append([strands[i][j], strands[i + 1][j]])
        if len(bonds) != num_bonds:
            bonds = []

        elif format_bonds(bonds)[-1][-1] == nbeads:
            print('index error: last index == nbeads')
            bonds = []

        else:
            return tuple((nbeads, len(bonds), str(format_bonds(bonds))))


def get_disulf_row(res_rows, num_bonds):
    """
    Function for generating  the disulfide bridge bonds

    Parameters
    ----------
    num_bonds: int. The number of bonds in this structure
    res_rows: set. Contains the structure row values (nbeads, nbonds, bead_index)

    Returns
    -------
    None. Adds new row values of nbeads, nbonds and bead_index for a disulfide structure to the rows set.
    """

    minbeads, maxbeads = num_bonds * 5, num_bonds * 5 + 5
    nbeads = ri(minbeads, maxbeads)
    possible_indices = set(el for el in range(nbeads))

    bonds = []
    while num_bonds:
        bond = sample(possible_indices, k=2)
        while abs(bond[0] - bond[1]) < 3:
            bond = sample(possible_indices, k=2)

        bonds.append(bond)
        possible_indices -= set(bond)
        num_bonds -= 1

    res_rows.add((nbeads, len(bonds), str(bonds)))


def get_mixed_row(num_bonds):
    """
    Function for getting mixed structure bonds
    Parameters
    ----------
    num_bonds: int. The number of bonds in this structure

    Returns
    -------
    tuple
    """

    minbeads, maxbeads = max(15, num_bonds * 2 + 2), min(40, num_bonds * 5)
    nbeads = ri(minbeads, maxbeads)
    bpb = [ri(1, 2) for _ in range(nbeads)]  # randomize the number of bonds possible for each bead
    possible_indices = [i for i, bp in enumerate(bpb) for _ in range(bp)]

    bonds = []
    while num_bonds and len(possible_indices) > 1:
        bond = sample(possible_indices, k=2)
        while abs(bond[0] - bond[1]) < 3:
            bond = sample(possible_indices, k=2)
        bond.sort()

        if bond not in bonds:
            bonds.append(bond)
            possible_indices.remove(bond[0])
            possible_indices.remove(bond[1])
            num_bonds -= 1

    if format_bonds(bonds)[-1][-1] == nbeads: print('index error: last index == nbeads')
    return tuple((nbeads, len(bonds), str(format_bonds(bonds))))


def get_rows(num_total=0, num_alpha=0, num_beta=0, num_disulf=0, num_mixed=0, num_bonds=1, all_rows=None):
    """
    Function takes in the number of structures required for each category and returns a set containing all of their
    indices, size and number of bonds

    Parameters
    ----------
    all_rows: the set of all rows ontaining the strcutural information
    num_bonds: int. The number of bonds in this structure
    num_alpha: int. The number of alpha helices to make
    num_beta: int. The number of beta sheets to make
    num_disulf: int. The number of disulfide bridge structures to make
    num_mixed: int. The mixed strcutures to generate

    Returns
    -------
    all_rows: set
    """

    # Chekc if all_rows has any elements already if not create a set
    if all_rows is None:
        all_rows = set()

    pre_len = len(all_rows)  # find the lenght of all_rows before strucutres put in it

    # Following code makes the different sets of strcutures
    alpha_rows = set()

    # it is possible to exhaust all unique alpha helix combinations so just try num_alpha * 100 times
    for _ in range(num_alpha * 100):
        alpha_rows.add(get_alpha_helix_row(num_bonds=num_bonds))
        if len(alpha_rows) == num_alpha:
            all_rows.update(alpha_rows)
            break
    all_rows.update(alpha_rows)

    beta_rows = set()

    # Can obtain all possible betasheet rows so use while loop
    while len(beta_rows) < num_beta:
        beta_rows.add(get_beta_row(num_bonds=num_bonds))
    all_rows.update(beta_rows)

    for _ in range(num_disulf):
        get_disulf_row(res_rows=all_rows, num_bonds=num_bonds)

    mixed_rows = set()
    while len(mixed_rows) < (num_mixed):
        mixed_rows.add(get_mixed_row(num_bonds=num_bonds))
    all_rows.update(mixed_rows)

    # In case some of the alpha helix strcutres quota could nto be filled pad with mixed structures
    while len(all_rows) < (num_total + pre_len):
        all_rows.add(get_mixed_row(num_bonds=num_bonds))

    print(len(alpha_rows), len(beta_rows), num_mixed, len(all_rows) - (len(alpha_rows) + len(beta_rows)), len(all_rows))
    return all_rows


def get_json(rows, tag=''):
    """
    Parameters
    ----------
    rows: set. Contains the structure row values (nbeads, nbonds, bead_index)
    
    Returns
    -------
    None. Creates the json file with all the structure information
    """
    # Get experiments
    nbeads, max_nbonds, nonlocal_bonds = zip(*rows)

    # Get some new column values for each row
    rc = [ri(15, 23) / 10 for _ in range(len(rows))]  # Get some random rc values for each bond
    nonlocal_bonds = [json.loads(el) for el in nonlocal_bonds]  # Change field values' type from str back to list.
    transient_bonds = nonlocal_bonds  # Set transient equal to non-local bonds
    config_out = [2 ** el for el in max_nbonds]  # Get config_out for each row from its max_nbonds value

    # Zip together the different columns, now including the new ones
    exp_vals = zip(rc, nonlocal_bonds, transient_bonds, config_out, nbeads, max_nbonds)

    # Create dictionary to store the new row values
    dict_for_json = [dict(zip(("rc", "nonlocal_bonds", "transient_bonds", "config_out", "nbeads", "max_nbonds"), row))
                     for row in exp_vals]

    # Put in the common values
    for el in dict_for_json:
        for key in common:
            el[key] = common[key]

    # Write each into json files
    for n, row in enumerate(dict_for_json):
        # Check what the last written
        # with open(f'../training_set_{max_nbonds[0]}_{tag}/train_{max_nbonds[0]}bonds_{n}.json', 'w') as jf:
        last_index = len([el for el in os.listdir() if (el.endswith('.json') and el.startswith('train_'))])
        with open(f'train_{max_nbonds[0]}bonds_{n + last_index}.json', 'w') as jf:
            json.dump(row, jf)

    # return dict_for_json


# Put rows in csv
# noinspection PyShadowingNames
def update_csv(csv_file, initialize_file, rows):
    """
    Parameters
    ----------
    rows: set. The set containing the row information for all the structures.
    csv_file: str. The csv file name to append to.
    initialize_file: int. If true then create csv file with required columns (cols).

    Returns
    -------
    None. Appends row values from numbeads, numbonds and bead_index compilations to the given csv file.
    """
    if initialize_file:
        with open(csv_file, 'w', newline='') as cf:
            csv_writer = csv.writer(cf, delimiter=',')
            csv_writer.writerow(cols)

    with open(csv_file, 'a', newline='') as cf:
        csv_writer = csv.writer(cf, delimiter=',')
        for row in rows:
            csv_writer.writerow(row)


def existing_rows():
    pre_rows = set()
    for el in os.listdir():
        if el.endswith('.json'):
            with open(el, 'r') as ij:
                data = json.load(ij)
            nb = data['nonlocal_bonds']
            nb = format(nb)
            pre_rows.add(tuple((data['nbeads'], len(nb), str(nb))))


def main(nbonds=3, num_rows=num_rows, tag='_v1'):
    get_json(
        get_rows(

            # Use dictionaru num_rows to get number of structures to generate for each bond type
            # For this number of bonds
            num_mixed=num_rows[nbonds]['mixed'],
            num_alpha=num_rows[nbonds]['alpha'],
            num_beta=num_rows[nbonds]['beta'],

            # nbonds and existing strucs (avoid repeats)
            num_bonds=nbonds,
            all_rows=existing_rows()
        ),
        # version number
        tag=tag
    )


if __name__ == "__main__":
    nbonds = int(sys.argv[1])
    tag = sys.argv[2]
    main(nbonds=nbonds, tag=tag)
