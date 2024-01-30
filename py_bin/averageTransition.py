
import glob
from pandas import read_csv
import csv
import numpy as np
def compute_averages():
    src_bias, src_mfpt = '*/diff_s_bias.csv', '*/mfpt.csv'

    diff_sbias = []
    n_ratio = []
    inner_mfpts = []
    outer_mfpts = []
    # average diff_s_bias
    for csv_file in glob.glob(src_bias):
        sbias_data = read_csv(csv_file, header=None)
        diff_sbias.append(sbias_data.iloc[0, 2])
        n_ratio.append(np.exp(sbias_data.iloc[0, 2]))

    for csv_file in glob.glob(src_mfpt):
        mfpts_data = read_csv(csv_file)
        inner_mfpts.append(mfpts_data.iloc[:, 3].values[0])
        outer_mfpts.append(mfpts_data.iloc[:, 4].values[0])

    output = []
    for data in [diff_sbias, n_ratio, inner_mfpts, outer_mfpts]:
        mean, var = np.mean(data), np.var(data)
        percent_rel_e = 100 * np.sqrt(var) / mean
        output.append([mean, var, percent_rel_e])


    csv_name = 'summary_data.csv'
    with open(csv_name, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)



def main():
    compute_averages()

if __name__ == '__main__':
    main()
