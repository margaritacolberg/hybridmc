#!/usr/bin/env python3
# Copyright (c) 2018-2021 Margarita Gladkikh
# SPDX-License-Identifier: BSD-3-Clause

import glob
import h5py
import json
import numpy as np
import os
from scipy import sparse

def main():
    h5_name = 'dataset.h5'
    h5_out = h5py.File(h5_name, 'w')

    src = '*/*/dataset.json'

    for json_in in glob.glob(src):
        if os.path.exists(json_in):
            print("Input json:", json_in)

            group_name = os.path.split(os.path.dirname(json_in))[0]

            with open(json_in, 'r') as input_json:
                data = json.load(input_json)

            nstates = len(data)

            f_map = []
            f_map_180 = []
            inputs = []
            outputs = [] 
            for i in range(nstates):
                nbeads = data[str(i)]['nbeads']
                nl_bonds = data[str(i)]['nonlocal_bonds']
                p_bonds = data[str(i)]['permanent_bonds']
                rc = data[str(i)]['rc']
                s_bias = data[str(i)]['s_bias']

                f_map_i, f_map_i_180, rc, nbonds, frac_bonds_on, frac_bonds_off, frac_beads_on, frac_beads_off = input_features(nbeads, nl_bonds, p_bonds, rc)

                f_map.append(f_map_i)
                f_map_180.append(f_map_i_180)

                inputs.append([rc, nbonds, frac_bonds_on, frac_bonds_off, frac_beads_on, frac_beads_off])
                outputs.append(s_bias)

            protein = h5_out.create_group(group_name)
            f_data = protein.create_dataset('f_map', data=f_map)
            f_180_data = protein.create_dataset('f_map_180', data=f_map_180)
            i_data = protein.create_dataset('inputs', data=inputs)
            o_data = protein.create_dataset('outputs', data=outputs)

def input_features(nbeads, nl_bonds, p_bonds, rc):
    nbonds = len(nl_bonds)
    nbonds_on = len(p_bonds)
    nbonds_off = nbonds - nbonds_on
    frac_bonds_on = nbonds_on / nbonds
    frac_bonds_off = nbonds_off / nbonds

    beads_in_p_bonds = 0
    for i in range(nbonds_on):
        beads_in_p_bonds = beads_in_p_bonds + p_bonds[i][1] - p_bonds[i][0]

    beads_in_nl_bonds = 0
    for i in range(nbonds):
        beads_in_nl_bonds = beads_in_nl_bonds + nl_bonds[i][1] - nl_bonds[i][0]

    frac_beads_on = beads_in_p_bonds / beads_in_nl_bonds
    frac_beads_off = 1.0 -  beads_in_p_bonds / beads_in_nl_bonds

    # nearest neighbors
    upper_1st_diag = np.ones(nbeads-1)
    lower_1st_diag = np.ones(nbeads-1)
    # next-nearest neighbors
    upper_2nd_diag = np.ones(nbeads-2)
    lower_2nd_diag = np.ones(nbeads-2)

    f_diags = [upper_1st_diag, lower_1st_diag, upper_2nd_diag, lower_2nd_diag]
    # feature map of local and nonlocal bonds
    f_map = sparse.diags(f_diags, [1, -1, 2, -2], dtype=int).toarray()

    if p_bonds:
        for k in range(len(p_bonds)):
            i = p_bonds[k][0]
            j = p_bonds[k][1]
            f_map[i][j] = 1
            f_map[j][i] = 1

    f_map_90 = np.rot90(f_map)
    f_map_180 = np.rot90(f_map, 2)
    f_map_270 = np.rot90(f_map, 3)

    return f_map, f_map_180, rc, nbonds, frac_bonds_on, frac_bonds_off, frac_beads_on, frac_beads_off

if __name__ == '__main__':
    main()
