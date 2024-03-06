import glob
from csv import DictWriter
import numpy as np
from pandas import read_csv
#from scipy.stats import bootstrap
#from sklearn.utils import resample as b_strap


def summarize_transition(D_mean=0.044062, D_std_perc=5):
    # define patterns for diff s bias and mfpt files
    src_bias, src_mfpt = '*/diff_s_bias.csv', '*/mfpt.csv'

    # intitialize lists to hold sbias and mfpts for the transition
    diff_sbias = []
    inner_mfpts = []
    outer_mfpts = []
    eps_vals = []
    pop_vals = []
     

    # average diff_s_bias compilation
    for csv_file in glob.glob(src_bias):
        sbias_data = read_csv(csv_file, header=None)
        diff_sbias.append(sbias_data.iloc[0, 2])
        eps_vals.append(3.0)
        pop_vals.append( np.exp( sbias_data.iloc[0,2]) )

    # mfpt compilation
    for csv_file in glob.glob(src_mfpt):
        mfpts_data = read_csv(csv_file)
        inner_mfpts.append(mfpts_data.iloc[:, 3].values[0])
        outer_mfpts.append(mfpts_data.iloc[:, 4].values[0])

    D_vals = np.random.normal(D_mean, D_std_perc/100.*D_mean, len(diff_sbias))
    std_D_test = np.std(D_vals)
    print(' test of std error of D is ', std_D_test/np.mean(D_vals) * 100, ' compare to ', D_std_perc)
    # put values together in a list of dictionaries
    output = [
        {'ID': 'sbias', 'vals': diff_sbias},
        {'ID': 'relative_population', 'vals': np.exp(diff_sbias)},
        {'ID': 'inner_mfpts', 'vals': inner_mfpts},
        {'ID': 'outer_mfpts', 'vals': outer_mfpts}
    ]

    # for each dict in output find relevant stats
    for data in output:
        data['count'], data['mean'], data['std'] = len(data['vals']), np.mean(data['vals']), np.std(data['vals'])
        data['percent_rel_e'] = 100 * data['std'] / data['mean']
        data['conf_int'] = f"{data['mean'] - 1.96 * data['std']}, {data['mean'] + 1.96 * data['std']}"
        data.pop('vals')

    # calculate array of predicted K values
    K_array = []
    for i in range(len(diff_sbias)):
        k_val = compute_rate_constant(D_vals[i],3.0, pop_vals[i], inner_mfpts[i], outer_mfpts[i])
        K_array.append(k_val)


    # obtain the rate constant (K) summary data
    K_mean, K_std, K_lower_conf, K_upper_conf = compute_K_summary(output, npoints=data['count'], D_mean=D_mean, D_std_perc=D_std_perc)

    K_mean2 = np.mean(K_array)
    K_std2 = np.std(K_array)
    print('direct average of K_array gives K_mean = ', K_mean2, ' and intervals (',
          K_mean2 - 1.96*K_std2, ',', K_mean2 + 1.96*K_std2, ')')

    # add the rate constant information to the output as well
    output.append({
        'ID': 'rate_constant',
        'count': output[0]['count'],
        'mean': K_mean, 'std': K_std,
        'percent_rel_e': 100 * K_std / K_mean,
        'conf_int': f'{K_lower_conf}, {K_upper_conf}'
    })

    # write the dict with stats to a csv file
    csv_name = 'summary_data.csv'
    with open(csv_name, 'w') as output_csv:
        writer = DictWriter(output_csv, fieldnames=output[0].keys())
        writer.writeheader()
        writer.writerows(output)

    print('Completed compiling and writing transition stats for sbias and mfpt')


def get_values(output, ID):
    for el in output:
        if el['ID'] == ID:
            return el['mean'], el['std']


def compute_rate_constant(D, eps, rel_population, inner_mfpt, outer_mfpt):
    return D * (1 + np.exp(-eps) * rel_population) / (outer_mfpt + np.exp(-eps) * rel_population * inner_mfpt)


def compute_K_summary(summary_data, npoints, D_mean=0.0403, D_std_perc=5, eps=3):
    rel_population, rel_population_std = get_values(summary_data, 'relative_population')
    inner_mfpt, inner_mfpt_std = get_values(summary_data, 'inner_mfpts')
    outer_mfpt, outer_mfpt_std = get_values(summary_data, 'outer_mfpts')
    #print('npoints is ',npoints)

    D_std = D_mean * D_std_perc / 100

    # rate constant K
    K = compute_rate_constant(D_mean, eps, rel_population, inner_mfpt, outer_mfpt)

    # canonical probability for given energy eps
    eps_population = np.exp(-eps) * rel_population

    # inner mfpt derivative prefactor
    alpha_inner = -eps_population / (outer_mfpt + inner_mfpt * eps_population)

    # outer mfpt derivative prefactor
    alpha_outer = -1 / (outer_mfpt + inner_mfpt * eps_population)

    # diffusion coefficient (D) derivative prefactor
    alpha_D = 1 / D_mean

    alpha_rel_population = (np.exp(-eps) / (1 + eps_population)) - \
                           np.exp(-eps) * inner_mfpt / (outer_mfpt + eps_population * inner_mfpt)

    rel_error_sq = 0
    for alpha, std in [(alpha_D, D_std), (alpha_rel_population, rel_population_std), (alpha_inner, inner_mfpt_std),
                       (alpha_outer, outer_mfpt_std)]:
        rel_error_sq += (alpha * std) ** 2

    # standard deviation of averages gives the standard error of estimated K
    K_std = np.sqrt(rel_error_sq) * K
    K_std_error = K_std

    K_lower_conf, K_upper_conf = K - 1.96 * K_std_error, K + 1.96 * K_std_error

    print(f'The rate constant is {K} with std {K_std}')

    print(f'The estimated 95% confidence intervals are {K_lower_conf}, {K_upper_conf}')

    return K, K_std, K_lower_conf, K_upper_conf


def main(args):
    summarize_transition(args.D_mean, args.D_std_perc)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--D_mean', type=float, default = 0.046, help='mean diffusion rate')
    parser.add_argument('--D_std_perc', type=float, default=5.0, help='percent error in the diffusion rate')

    args = parser.parse_args()

    main(args)
