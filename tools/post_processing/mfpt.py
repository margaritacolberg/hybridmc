#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# mfpt.py determines the mean first passage times using spline interpolation
# and their errors using the bootstrap method; the outer first passage time
# (fpt) corresponds to the first time at which two beads at a distance r > rc
# diffuse to rc, and the inner fpt corresponds to beads located at r < rc
#
# example of how to run:
# python ../tools/mfpt.py hybridmc_0_000_001.json hybridmc_0_000_001.h5 10 hybridmc_0_000_001.csv --layers=True
#
# note that in the above example, hybridmc_0_000_001.json is the json file for
# the last step of the staircase (ie. where the width of the step is between rh
# and rc)

from . import minimize, matrix_element
import multiprocessing
import csv
import h5py
import json
import numpy as np
import nlopt
from scipy import integrate
from scipy import interpolate
from scipy import sparse
from sklearn import utils

import os
import re
import glob

def process_mfpts():
    output = []
    for file_path in os.listdir('.'):
        if re.search('hybridmc_[0-9]+_([01]+)_([01]+).csv', file_path):
            bits_i, bits_j, t_ind, int_i, int_j = matrix_element.get_state_data(file_path)
            inner_fpt, outer_fpt = matrix_element.get_fpt(file_path, t_ind)
            layer = str(bits_i).count('1')
            output.append([bits_i, bits_j, layer, inner_fpt, outer_fpt])

    output.sort(key=lambda x: x[2])
    csv_out = 'mfpt.csv'

    output.insert(0, ['state i bits', 'state j bits', 'layer', 'inner fpt',
                      'outer fpt'])

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def get_mfpt():
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
        pool.starmap(fpt_write, names)


def fpt_write(json_in, hdf5_in, nboot, csv_out, layers):
    try:
        fpt_write_wrap(json_in, hdf5_in, nboot, csv_out, layers)
    except:
        print(json_in, 'FAILED')
        raise


def fpt_write_wrap(json_in, hdf5_in, nboot, csv_out, layers):
    print('input json:', json_in)
    with open(json_in, 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']
    beta = 1.0
    nknots = 12

    t_bonds = data['transient_bonds']
    p_bonds = data['permanent_bonds']
    nl_bonds = data['nonlocal_bonds']

    rc_transient = t_bonds[-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    if p_bonds:
        p_ind = [i for i, bond in enumerate(nl_bonds) for j in
                 range(len(p_bonds)) if bond == p_bonds[j]]
    else:
        p_ind = []

    with h5py.File(hdf5_in, 'r') as f:
        dist = f['dist'][:]

    print('number of distances =', len(dist))

    dist = dist.T
    nbonds = len(nl_bonds)

    # if running layers for entire folding process, run bootstrap samples for
    # a specific transition to determine variance
    if layers and nboot > 0:
        if t_ind == 1 and len(p_ind) == 1 and p_ind[0] == 0:
            run_bootstrap = True
        else:
            run_bootstrap = False
    # if running one transition only (ie. for testing),
    else:
        if not layers and nboot > 0:
            run_bootstrap = True
        else:
            run_bootstrap = False

    output = []
    t_on = []
    t_off = []
    for j in range(len(dist[t_ind])):

        if dist[t_ind][j] < rc_transient:
            t_on.append(dist[t_ind][j])

        else:
            t_off.append(dist[t_ind][j])

        t_on = np.array(t_on)
        t_off = np.array(t_off)

        # inner and outer fpt for transient bond turned on
        fpt_on = fpt_per_bead_pair(t_on, nknots, beta, rh, rc_transient, True)[0]
        fpt_off = fpt_per_bead_pair(t_off, nknots, beta, rc_transient, np.max(t_off), False)[0]

        if run_bootstrap:
            fpt_on_var = fpt_var(t_on, nknots, beta, rh, rc_transient, True, nboot)
            fpt_off_var = fpt_var(t_off, nknots, beta, rc_transient, np.max(t_off), False, nboot)
            output += [fpt_on_var, fpt_off_var]

        output = [t_ind, fpt_on, fpt_off]

    with open(csv_out, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def A_matrix(x):
    upper_A_diag = np.zeros(x.size - 1)
    A_diag = np.ones(x.size)
    lower_A_diag = np.zeros(x.size - 1)

    h_lower = x[1:-1] - x[:-2]
    h_upper = x[2:] - x[1:-1]

    lower_A_diag[:-1] = h_lower
    A_diag[1:-1] = 2 * (h_upper + h_lower)
    upper_A_diag[1:] = h_upper

    A_diags = [A_diag, lower_A_diag, upper_A_diag]
    A = sparse.diags(A_diags, [0, -1, 1]).toarray()

    return np.linalg.inv(A)


def spline(x_knot, y_knot, x_input):
    spl = interpolate.CubicSpline(x_knot, y_knot, bc_type='natural')
    coeffs = np.flip(spl.c, 0)
    a, b, c, d = coeffs[0], coeffs[1], coeffs[2], coeffs[3]

    # list whose elements are the index of each knot to the left of each x
    # input
    ind = np.searchsorted(x_knot, x_input, side='left')

    if ind[0] == 0:
        ind[1:] -= 1
    else:
        ind -= 1

    f = []
    for i in range(len(ind)):
        del_x = x_input[i] - x_knot[ind[i]]
        coeff_i = np.column_stack((a[ind[i]], b[ind[i]], c[ind[i]], d[ind[i]]))
        f.append((coeff_i[:, 0] + coeff_i[:, 1] * del_x + coeff_i[:, 2] * del_x ** 2
                  + coeff_i[:, 3] * del_x ** 3)[0])

    return np.array(f)


# norm factor
def integrate_spline(x_knot, y_knot):
    output = []
    f = lambda x: np.exp(-spline(x_knot, y_knot, x))

    for i in range(len(x_knot) - 1):
        output.append(integrate.quadrature(f, x_knot[i], x_knot[i + 1], tol=1e-6,
                                           rtol=1e-6)[0])

    return sum(output)


def fun(x, y, x_i):
    integral = integrate_spline(x, y)

    f = np.sum(spline(x, y, x_i))
    count = len(x_i)

    if count > 0:
        f /= (count * 1.0)

    f += np.log(integral)

    return f


def dB_matrix(x):
    upper_dB_diag = np.zeros(x.size - 1)
    dB_diag = np.zeros(x.size)
    lower_dB_diag = np.zeros(x.size - 1)

    h_lower = 3 / (x[1:-1] - x[:-2])
    h_upper = 3 / (x[2:] - x[1:-1])

    lower_dB_diag[:-1] = h_lower
    dB_diag[1:-1] = -h_upper - h_lower
    upper_dB_diag[1:] = h_upper

    dB_diags = [dB_diag, lower_dB_diag, upper_dB_diag]

    return sparse.diags(dB_diags, [0, -1, 1]).toarray()


def dspline_coefficients(x, i):
    A = A_matrix(x)
    dB_i = dB_matrix(x)[:, i]

    # solve the linear system Adc = dB to get vector dc whose entries are the
    # dc coefficients of the spline
    dc_i = np.dot(A, dB_i)
    h = x[1:] - x[:-1]
    db = (-h / 3) * (2 * dc_i[:-1] + dc_i[1:])

    if i >= 1:
        db[i - 1] += 1 / h[i - 1]

    if i < len(db):
        db[i] -= 1 / h[i]

    dd = (dc_i[1:] - dc_i[:-1]) / (3 * h)

    dc_i = dc_i[:-1]

    return db, dc_i, dd


def dspline(x_knot, x_input, ind):
    db, dc, dd = dspline_coefficients(x_knot, ind)

    i = np.searchsorted(x_knot, x_input, side='left')

    if i[0] == 0:
        i[1:] -= 1
    else:
        i -= 1

    del_x = x_input - x_knot[i]
    coeff_mat = np.column_stack((db[i], dc[i], dd[i]))
    f = coeff_mat[:, 0] * del_x + coeff_mat[:, 1] * del_x ** 2 + coeff_mat[:,
                                                                 2] * del_x ** 3

    # add da / dy_i to spline, if it exists
    f[i == ind] += 1.0

    return np.array(f)


def integrate_dspline(x_knot, y_knot, ind):
    output = []
    f = lambda x: -dspline(x_knot, x, ind) * np.exp(-spline(x_knot, y_knot, x))

    for i in range(len(x_knot) - 1):
        output.append(integrate.quadrature(f, x_knot[i], x_knot[i + 1], tol=1e-6,
                                           rtol=1e-6)[0])

    return sum(output)


def dfun(x, y, x_i):
    f_integral = integrate_spline(x, y)

    n_derivs = len(y)
    count = len(x_i)

    df = []
    for i in range(n_derivs):
        df_integral = integrate_dspline(x, y, i)
        df_i = np.sum(dspline(x, x_i, i))

        if count > 0:
            df_i /= (count * 1.0)

        df_i += (1 / f_integral) * df_integral
        df.append(df_i)

    return np.array(df)


def minimize_f(x, x_i, y0):
    f = lambda y: fun(x, y, x_i)
    df = lambda y: dfun(x, y, x_i)

    # min_f2 = optimize.minimize(f, x0=y0, method='BFGS', jac=df,
    # options={'gtol': 0.000001, 'disp': True})
    min_f = minimize.minimize(f, y0, method=nlopt.LD_LBFGS, jac=df, ub=None,
                              lb=None)

    return min_f


def cdf_at_x(x_knot, y_knot, x, norm):
    output = []
    f = lambda t: np.exp(-spline(x_knot, y_knot, t))

    for i in range(len(x_knot) - 1):
        if (x_knot[i + 1] < x) and (x_knot[i] < x):
            output.append(integrate.quadrature(f, x_knot[i], x_knot[i + 1],
                                               tol=1e-6, rtol=1e-6)[0])
        elif (x_knot[i + 1] > x) and (x_knot[i] < x):
            output.append(integrate.quadrature(f, x_knot[i], x, tol=1e-6,
                                               rtol=1e-6)[0])

    return sum(output) / norm


def fpt_integrand(x_knot, y_knot, x, beta, state):
    if type(state) != bool:
        raise TypeError('the variable state must be bool')

    norm = integrate_spline(x_knot, y_knot)
    f = spline(x_knot, y_knot, x)

    pdf = np.exp(-beta * f) / norm

    cdf = []
    for i in range(len(x)):
        cdf.append(cdf_at_x(x_knot, y_knot, x[i], norm))

    cdf = np.array(cdf)

    integrand = 0.0
    if state is True:
        # fpt is inner
        integrand = cdf * cdf / pdf
    else:
        integrand = (1.0 - cdf) * (1.0 - cdf) / pdf

    return integrand


def fpt(x_knot, y_knot, beta, xmin, xmax, state):
    f = lambda x: fpt_integrand(x_knot, y_knot, x, beta, state)

    return integrate.quadrature(f, xmin, xmax, tol=1e-6, rtol=1e-6,
                                maxiter=100)


def fpt_per_bead_pair(dist_vec, nknots, beta, min_dist, max_dist, state):
    dist_vec.sort()

    x = np.linspace(min_dist, max_dist, nknots)
    # initial guess for y-coordinates of knots
    y = np.zeros(nknots)

    y = minimize_f(x, dist_vec, y)

    return fpt(x, y, beta, min_dist, max_dist, state)


def fpt_var(dist_vec, nknots, beta, min_dist, max_dist, state, nboot):
    sample_size = int(len(dist_vec) / 2)

    boot_fpt = []
    for i in range(nboot):
        print('bootstrap round', i)
        boot_fpt_i = utils.resample(dist_vec, n_samples=sample_size,
                                    random_state=i)
        boot_fpt.append(fpt_per_bead_pair(boot_fpt_i, nknots, beta, min_dist,
                                          max_dist, state)[0])

    return np.var(boot_fpt, ddof=1)


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))