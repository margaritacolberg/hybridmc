import csv
import json
import glob
import os


class input_features:

    def __init__(self, rc, nbeads, nl_bonds, p_bonds):
        self.rc = rc
        self.nbeads = nbeads
        self.nl_bonds = nl_bonds
        self.p_bonds = p_bonds

        self.nbonds = len(nl_bonds)
        self.nbonds_on = len(p_bonds)
        self.frac_bonds_on = self.nbonds_on / self.nbonds
        
        self.beadind_p_diff = self.accum_diff(self.p_bonds)
        self.beadind_nl_diff = self.accum_diff(self.nl_bonds)

        self.beadind_p_invar_sums = self.accum_invar_sum(self.p_bonds)
        self.beadind_nl_invar_sums = self.accum_invar_sum(self.nl_bonds)

        self.beadind_p_sums = self.accum_sum(self.p_bonds)
        self.beadind_nl_sums = self.accum_sum(self.nl_bonds)

    def accum_diff(self, bondlist):
        return sum([el[1] - el[0] for el in bondlist]) / self.nbeads

    def accum_invar_sum(self, bondlist):
        return sum([abs(
                (el[1] + el[0] + 1) / self.nbeads - 1
        ) for el in bondlist])

    def accum_sum(self, bondlist):
        return sum([(el[1] + el[0]) / self.nbeads for el in bondlist])

    @staticmethod
    def feature_names():
        return [
            'rc', 'nbonds', 'nbeads', 'frac_bonds_on',
            'beadind_p_diff', 'beadind_nl_diff',
            'beadind_p_sums', 'beadind_nl_sums'
            #'beadind_p_invar_sums', 'beadind_nl_invar_sums'
        ]

    def features(self):
        return [
            self.rc, self.nbonds, self.nbeads, self.frac_bonds_on,
            self.beadind_p_diff, self.beadind_nl_diff,
            #self.beadind_p_sums, self.beadind_nl_sums,
            self.beadind_p_invar_sums, self.beadind_nl_invar_sums
        ]


def write_csvrows_fromjson(json_in, cw):
    with open(json_in, 'r') as input_json:
        data = json.load(input_json)

    ref = 0
    if '0' in data:
        ref = data['0']['s_bias']
        del data['0']

    for i in data:
        row = input_features(
            nbeads=data[i]['nbeads'],
            nl_bonds=data[i]['nonlocal_bonds'],
            p_bonds=data[i]['permanent_bonds'],
            rc=data[i]['rc']
        ).features()

        row.extend([data[i]['s_bias'] - ref, data[i]['permanent_bonds'], data[i]['nonlocal_bonds'], os.path.dirname(json_in)])
        cw.writerow(row)


def make_csv(file, path, fieldnames, mode='w'):

    with open(file, mode, newline='') as cf:
        csv_writer = csv.writer(cf, delimiter=',')

        if mode == 'w':
            csv_writer.writerow(fieldnames)

        for json_in in glob.glob(path):
            print("Input json:", json_in)
            write_csvrows_fromjson(json_in=json_in, cw=csv_writer)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-cf', '--csv_file', type=str, default='9bond_data.csv')
    parser.add_argument('-s', '--src', type=str, default='')
    parser.add_argument('-m', '--mode', type=str, default='w')
    args = parser.parse_args()

    make_csv(
        args.csv_file, args.src, input_features.feature_names() + ['s_bias', 'p_bonds', 'nl_bonds', 'structure_id'],
        mode=args.mode)

