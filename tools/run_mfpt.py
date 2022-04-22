#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

import mfpt

import os
import multiprocessing
import glob


src = '*.json'

names = []
for file_path in glob.glob(src):
    file_name = os.path.basename(file_path)
    name = os.path.splitext(file_name)[0]

    json_in = name + '.json'
    hdf5_in = name + '.h5'
    csv_out = name + '.csv'

    if os.path.exists(csv_out):
        continue

    names.append((json_in, hdf5_in, 0, csv_out, False))

with multiprocessing.Pool() as pool:
    pool.starmap(mfpt.fpt_write, names)
