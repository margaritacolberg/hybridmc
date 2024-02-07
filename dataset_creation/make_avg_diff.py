import csv
import h5py
import os
from numpy import mean


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))


# get some information like nbonds from the directory name; might use some in csv name scheme
def get_nameinfo():
    """
    The function makes use of the directory name to extract some information about the simulation parameters
    Returns
    -------
    Tupe containing the nbonds, the tag, array number for the job and the path for the json file in here describind the
    simulated structure (jpath)
    """
    cp = os.getcwd()
    dir_split = os.path.splitext(cp)[0].rsplit(sep='_')
    nbonds, tag, array_num = dir_split[-3], dir_split[-2], dir_split[-1]
    jpath = f"{os.path.dirname(cp)}/train_{nbonds}bonds_{array_num}.json"

    return int(nbonds), tag, array_num, jpath


def update_csv(csv_file, rows, cols=None):
    """
    Create a csv file where each row is described by the elements of rows
    Parameters
    ----------
    cols: List[str]. Default is None, header only written to csv file if provided.
    csv_file: str
    rows: List[Lists]

    Returns
    -------
    None
    """
    with open(csv_file, 'w', newline='') as cf:
        csv_writer = csv.writer(cf, delimiter=',')

        if cols is not None:
            csv_writer.writerow(cols)
            
        csv_writer.writerows(rows)


def get_dsb_rows_v2(patt, sep):
    output = []
    for path in os.listdir():
        if path.endswith(patt):
            config_in = os.path.splitext(path)[0].rsplit(sep=sep)[-2]
            config_out = os.path.splitext(path)[0].rsplit(sep=sep)[-1]
            with h5py.File(path, 'r') as f:
                output.append((config_in, config_out, float(f['s_bias'][0])))

    return output


def get_asb_rows(dsb_rows, nbonds):
    diff_s_bias = {(row[0], row[1]): row[2] for row in dsb_rows}

    # initially fully bonded
    work_list = [[True] * nbonds]
    bonded_config = format_bits(work_list[0])

    s_bias = dict()
    s_bias[bonded_config] = 0

    avg_s_bias = dict()
    while work_list:
        bits_in = work_list.pop(0)
        config_in = format_bits(bits_in)
        mean_s_bias_in = mean(s_bias[config_in])
        avg_s_bias[config_in] = mean_s_bias_in

        for i in range(nbonds):
            # if two beads are not bonded,
            if not bits_in[i]:
                # skip; do not form the bond
                continue

            # copy of initial state for every bonded pair of beads
            bits_out = bits_in.copy()
            # flip bit to break bond
            bits_out[i] = False

            config_out = format_bits(bits_out)
            diff_s_bias_out = diff_s_bias[(config_out, config_in)]

            s_bias_out = s_bias.setdefault(config_out, [])
            s_bias_out.append(diff_s_bias_out + mean_s_bias_in)

            # check if bond pattern already exists inside worklist;
            # only append unique bond patterns to worklist
            if not bits_out in work_list:
                work_list.append(bits_out)

    return sorted(avg_s_bias.items())


def main():
    
    # make diff_s_bias csv file
    dsb_rows = get_dsb_rows_v2(patt='.h5', sep='_')  # diff_s_bias rows
    update_csv(csv_file='diff_s_bias.csv',
               rows=dsb_rows)

    # make avg_s_bias csv file
    update_csv(csv_file='avg_s_bias.csv',
               rows=get_asb_rows(
                   dsb_rows=dsb_rows,
                   nbonds=int(os.path.splitext(os.getcwd())[0].rsplit(sep='_')[-3])))


if __name__ == '__main__':
    main()