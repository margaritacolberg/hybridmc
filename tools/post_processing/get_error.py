import glob
from csv import DictWriter, writer

import numpy as np
from pandas import read_csv


def summarize_transition():
    # define patterns for diff s bias and mfpt files
    src_bias, src_mfpt = '*/diff_s_bias.csv', '*/mfpt.csv'

    # intitialize lists to hold sbias and mfpts for the transition
    diff_sbias = []
    inner_mfpts = []
    outer_mfpts = []

    # average diff_s_bias compilation
    for csv_file in glob.glob(src_bias):
        sbias_data = read_csv(csv_file, header=None)
        diff_sbias.append(sbias_data.iloc[0, 2])

    # mfpt compilation
    for csv_file in glob.glob(src_mfpt):
        mfpts_data = read_csv(csv_file)
        inner_mfpts.append(mfpts_data.iloc[:, 3].values[0])
        outer_mfpts.append(mfpts_data.iloc[:, 4].values[0])

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

    # obtain the rate constant (K) summary data
    K_mean, K_std, K_lower_conf, K_upper_conf = compute_K_summary(output)

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


def compute_K_summary(summary_data, D_mean=0.046, D_std=0.0023, eps=3):
    rel_population, rel_population_std = get_values(summary_data, 'relative_population')
    inner_mfpt, inner_mfpt_std = get_values(summary_data, 'inner_mfpts')
    outer_mfpt, outer_mfpt_std = get_values(summary_data, 'outer_mfpts')

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

    K_std = np.sqrt(rel_error_sq) * K
    K_lower_conf, K_upper_conf = K - 1.96 * K_std, K + 1.96 * K_std

    print(f'The rate constant is {K} with std {K_std}')

    print(f'The 95% confidence intervals are {K_lower_conf}, {K_upper_conf}')

    return K, K_std, K_lower_conf, K_upper_conf


def compute_rate_constant(D, eps, rel_population, inner_mfpt, outer_mfpt):
    return D * (1 + np.exp(-eps) * rel_population) / (outer_mfpt + np.exp(-eps) * rel_population * inner_mfpt)


def main():
    summarize_transition()


if __name__ == '__main__':
    main()
