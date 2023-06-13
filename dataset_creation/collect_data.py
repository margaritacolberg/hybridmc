#!/usr/bin/env python3
# Copyright (c) 2018-2021 Margarita Gladkikh
# SPDX-License-Identifier: BSD-3-Clause

import csv
import glob
import json
import os
import re

def main():
    csv_in = 'avg_s_bias.csv'
    json_out = 'dataset.json'

    s_bias = []
    with open(csv_in, 'r') as input_csv:
        s_bias_data = csv.reader(input_csv, delimiter=',')

        for row in s_bias_data:
            s_bias.append(float(row[1]))

    src = 'hybridmc_*.json'

    config = []
    output = {}
    for json_in in glob.glob(src):
        if os.path.exists(json_in):
            print("Input json:", json_in)
            bits_from_json = re.search(r"_([01]+)_([01]+)\.json$", json_in)
            bits_in, bits_out = bits_from_json.groups()

            int_in = int(bits_in, base=2)
            int_out = int(bits_out, base=2)

            with open(json_in, 'r') as input_json:
                data = json.load(input_json)

            if int_in not in config:
                output[str(int_in)] = {
                    'nbeads': data['nbeads'],
                    'nonlocal_bonds': data['nonlocal_bonds'],
                    'permanent_bonds': data['permanent_bonds'],
                    'rc': data['rc'],
                    's_bias': s_bias[int_in],
                    }

                config.append(int_in)

            nbits = len(bits_in)
            if int_out == (2**nbits - 1) and int_out not in config:
                output[str(int_out)] = {
                    'nbeads': data['nbeads'],
                    'nonlocal_bonds': data['nonlocal_bonds'],
                    'permanent_bonds': data['nonlocal_bonds'],
                    'rc': data['rc'],
                    's_bias': s_bias[int_out],
                    }

                config.append(int_out)

    with open(json_out, 'w') as output_json:
        json.dump(output, output_json)

if __name__ == '__main__':
    main()
