#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause


import matrix_element
import csv
import os
import re


def if_stair(file_path, files):
    """
    Funtion to check if the file_path ahs associated staircased steps simulations results. If yes, then the paths of these
    intermediate steps' csv files with their mfpt information are compiled into a list and returned.

    :param file_path: the final step mfpt csv file path
    :param files: the list of files in this directory
    :return: list containing all the intermediate staircase mfpt csv files
    """

    # initialize output list merging all intermediate stair mfpt csv paths
    merge_paths = []
    # get the reference simulation id tag
    ref_sim_id = file_path.rstrip('.csv')

    # loop through csv files in directory
    for file in files:
        if file.endswith('.csv'):
            # obtain simulation tag for this file
            sim_id = file.rstrip('.csv')
            # if the reference tag and this are the same then add the filepath to the output list
            if ref_sim_id.split('_') == sim_id.split('_')[:-1]:
                merge_paths.append(file)

    return merge_paths


def compile_outer_fpt(stair_paths, t_ind):
    """
    Function to compile all the staircase outer fpts to one number that is returned
    :param stair_paths: list containing all the intermediate staircase mfpt csv files
    :param t_ind: the transient bond index for this simulation
    :return: the compiled outer_fpt value
    """
    outer_fpt = 0
    for file_path in stair_paths:
        outer_fpt += matrix_element.get_fpt(file_path, t_ind)[1]

    return outer_fpt


def compile_mfpt():
    """
    Function to compile all the mfpt csv files in the current directory into one csv file

    Returns
    -------
    None, but writes a csv file with the compiled mfpt data

    """
    output = []
    files = os.listdir()
    for file_path in files:

        if re.search('hybridmc_[0-9]+_([01]+)_([01]+).csv', file_path):

            bits_i, bits_j, t_ind, _, _ = matrix_element.get_state_data(file_path)
            inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)
            layer = file_path.split('_')[1]

            stair_paths = if_stair(file_path, files)
            if stair_paths:
                outer_fpt += compile_outer_fpt(stair_paths, t_ind)

            output.append([bits_i, bits_j, layer, inner_fpt, outer_fpt])

    output.sort(key=lambda x: x[2])
    csv_out = 'mfpt.csv'

    output.insert(0, ['state i bits', 'state j bits', 'layer', 'inner fpt',
                      'outer fpt'])

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


if __name__ == '__main__':
    compile_mfpt()
