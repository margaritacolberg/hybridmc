import csv
import json
import os

import h5py
import nlopt
import numpy as np
import copy
from scipy import sparse, interpolate, integrate
from scipy import stats
from scipy.special import kolmogorov
from sklearn import utils

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from .data_processing_helpers import set_defaults, if_stair

from ..post_processing import matrix_element

# tol, rtol, maxiter = 1e-3, 1e-3, 100
plot_data = False
write_data = False
BFGS = True

Adaptive = False
Adaptive_KS = False
q_cut = 0.5
Verbose_Convergence = True
Verbose_Bootstrap = True

#best_val = np.inf
#best_args = None

import numpy as np

#Define a fallback optimizer in case of exception
def fallback_optimizer(f, jac, x_init, lb, ub):
	opt = nlopt.opt(nlopt.LN_BOBYQA, len(x_init))
	obj_f = make_nlopt_f(f, jac)
	opt.set_min_objective(obj_f)
	if ub is not None:
		opt.set_upper_bounds(ub)
	if lb is not None:
		opt.set_lower_bounds(lb)

	opt.set_xtol_rel(1e-4)
	opt_results = opt.optimize(x_init)
	return opt_results

def minimize(f, x0, method, jac, ub, lb):
    #global best_args
    #global best_val

    dim = len(x0)

    try: 
        opt = nlopt.opt(method, dim)
        obj_f = make_nlopt_f(f, jac)
        opt.set_min_objective(obj_f)

        if ub is not None:
            opt.set_upper_bounds(ub)

        if lb is not None:
            opt.set_lower_bounds(lb)

        # stopping criteria
        #
        # function tolerance
        opt.set_ftol_rel(5e-5)
        # y-coordinate tolerance
        opt.set_xtol_rel(5e-5)



        opt_results = opt.optimize(x0)

    except nlopt.RoundoffLimited as e:
        print("NLopt failed due to roundoff errors. Switching to fallback optimizer.")
        x0 = np.zeros(dim) 
        opt_results = fallback_optimizer(f, jac, x0, lb, ub)
        print("Found minimum using fallback optimizer:", opt_results)  

    except Exception as e:
        print("NLopt failed due to an unexpected exception:", e)
        print("Switching to fallback optimizer.")
        x0 = np.zeros(dim)
        opt_results = fallback_optimizer(f, jac, x0, lb, ub)
        print("Found minimum using fallback optimizer:", opt_results)

    return opt_results


def make_nlopt_f(f, jac):
    def nlopt_f(x, grad):
        if grad.size > 0:
            grad[:] = jac(x)

        return f(x)

    return nlopt_f


def fpt_write(name):
    # case 1: if name is for an intermediate simulation do nothing
    if len(name.split('_')) == 5:
        #        print(f"Skipping {name}: it is an intermediate simulation")
        return

    else:
        # if name is for a final simulation, check if there are intermediate simulations
        stair_rc_list = if_stair(name, os.listdir())

        # case 2: if there are intermediate simulations, run the staircase mfpt
        if len(stair_rc_list) > 0:
            fpt_write_wrap_stair(name, stair_rc_list)

        # case 3: if there are no intermediate simulations, run the default mfpt
        else:
            fpt_write_wrap_default(name)


def get_dists(name, t_ind, max_length=25000):
    """
    Function to load the distances for the transient bond active now
    Parameters
    ----------
    name: str: the name of the simulation
    t_ind: int: the index of the transient bond active now
    max_length: int: the maximum number of distances to load

    Returns
    -------
    np.array: the distances for the transient bond active now
    """
    # load the distances for this simulation
    with h5py.File(f"{name}.h5", 'r') as f:
        # load only at most 25000 distances
        max_dist = min([max_length, len(f['dist'])])
        # initially the shape is (25000, nbonds). each column with dists for a transient bond
        dist = f['dist'][:max_dist]
    # get the transpose to make each row the dists for a transient bond
    dist = dist.T
    # pick out the distances for the transient bond active now
    dist_t_active = dist[t_ind]
    return dist_t_active


def fpt_write_wrap_default(name):
    print('input json:', f"{name}.json")
    with open(f"{name}.json", 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']

    # set some more default values for parameters in case they are not provided in json file
    set_defaults(data,
                 defaults={
                     "WL_sbias": 6.0,
                     "fail_max": 5,
                     "req_dists": 25000,
                     "rc_target_min_percentile": 0.10,
                     "nsteps_max": 1000,
                     "nboot": 0

                 })

    t_bonds = data['transient_bonds']
    p_bonds = data['permanent_bonds']
    nl_bonds = data['nonlocal_bonds']
    max_d = data['req_dists']
    nboot = data['nboot']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    if p_bonds:
        p_ind = [i for i, bond in enumerate(nl_bonds) for j in
                 range(len(p_bonds)) if bond == p_bonds[j]]
    else:
        p_ind = []

    # if running layers for entire folding process, run bootstrap samples for
    # a specific transition to determine variance
    if nboot:
        run_bootstrap = True
    else:
        run_bootstrap = False

    output = []

    # generate for non staircase simulations, the values of the distances within and outside the rc
    t_on = []
    t_off = []

    dist_t_active = get_dists(name, t_ind, max_d)
    for j in range(len(dist_t_active)):

        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])
        else:
            t_off.append(dist_t_active[j])

    t_on = np.array(t_on)
    t_off = np.array(t_off)

    # inner fpt for transient bond turned on
    xb_on, yb_on, fpt_on = fpt_per_bead_pair(t_on, rh, rc_transient, True, name)
    xb_off, yb_off, fpt_off = fpt_per_bead_pair(t_off, rc_transient, np.max(t_off), False, name)

    output_i = [t_ind, fpt_on, fpt_off]

    if run_bootstrap:
        fpt_on_std = fpt_std(t_on, rh, rc_transient, True, nboot, xb_on, yb_on)
        fpt_off_std = fpt_std(t_off, rc_transient, np.max(t_off), False, nboot, xb_off, yb_off)
        output_i += [fpt_on_std, fpt_off_std]
    else:
        output_i += [0.0, 0.0]

    output.append(output_i)

    with open(f"{name}.csv", 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def get_fpt_boot_distances(stair_output_name, boot_index):
    with open(f"{stair_output_name}.json", 'r') as input_json:
        data = json.load(input_json)

    t_bonds = data['transient_bonds']
    nl_bonds = data['nonlocal_bonds']
    max_d = data['req_dists']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    dist_t_active = get_dists(stair_output_name, t_ind, max_d)

    t_on = []
    t_off = []

    for j in range(len(dist_t_active)):
        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])
        else:
            t_off.append(dist_t_active[j])

    t_off = np.array(t_off)
    boot_fpt_i = utils.resample(t_off, random_state=boot_index)
    return boot_fpt_i


def fpt_outer_std(name, stair_rc_list, xknots_iteration, yknots_iteration, nboot):
    outer_fpt_array = []
    for boot in range(nboot):
        xknots_total = []
        yknots_total = []
        for j in range(len(stair_rc_list)):
            # So long as the rc is not at the target we append the intermediate rc to the output name
            if stair_rc_list[j] != stair_rc_list[-1]:
                stair_output_name = f'{name}_{stair_rc_list[j]}'
            else:
                stair_output_name = name

            if j == 0:
                rc_max = -1.
                # otherwise, the maximum rc is the rc of the previous step
            else:
                rc_max = stair_rc_list[j - 1]
                # if this is the first step, then the knots are the total knots that subsequent knots are added

            dist_vec_j = get_fpt_boot_distances(stair_output_name, boot_index=boot)

            xknots = xknots_iteration[j].copy()
            yknots = minimize_f(xknots_iteration[j], dist_vec_j, yknots_iteration[j])
            yknots -= yknots[-1]

            if j == 0:
                xknots_total = xknots.copy()
                yknots_total = yknots.copy()

                # otherwise, append the new knots to the total knots to build a smooth spline across the entire rc space
            else:
                # In the case of the x knots, simply append the new x knots to the total x knots less the first knot
                # of the total x knots list (this is repeated exactly in the new x knots list as its last value)
                xknots_total = np.append(xknots, xknots_total[1:])

                # In the case of the y_knot values:
                # 1. Do the same as the x knots AND
                # 2. Shift the new y knots by the first y knot value in the total knots list to avoid a discontinuity
                # in the spline fit at the junction of the two staircase steps
                yknots_total = np.append(yknots + yknots_total[0], yknots_total[1:])
        fpt_boot_val = fpt(xknots_total, yknots_total, xknots_total[0], xknots_total[-1], False)[0]

        if plot_data or write_data:
            plot_outer_integrand(xknots_total, yknots_total, f'{name}.{boot}')

        if Verbose_Bootstrap:
            print(' outer bootstrap ', boot, ' has value ', fpt_boot_val, ' for transition ', name)

        outer_fpt_array.append(fpt_boot_val)

    return np.std(outer_fpt_array)


def fpt_outer_stair(name, stair_rc_list):
    # get knots for the whole rc space to make one smooth spline fit across the outer range encapsulating all the rcs
    xknots_total = []
    yknots_total = []

    xknots_iteration = []
    yknots_iteration = []
    for j in range(len(stair_rc_list)):
        # So long as the rc is not at the target we append the intermediate rc to the output name
        if stair_rc_list[j] != stair_rc_list[-1]:
            stair_output_name = f'{name}_{stair_rc_list[j]}'
        else:
            stair_output_name = name

        # if the rc is the outermost rc (step 1), then the maximum rc can be the largest distance in the data file
        # indicate this by using -1 as the rc_max
        if j == 0:
            rc_max = -1.
        # otherwise, the maximum rc is the rc of the previous step
        else:
            rc_max = stair_rc_list[j - 1]
            # if this is the first step, then the knots are the total knots that subsequent knots are added

        nknots, xknots, yknots = knots_wrap(stair_output_name, rc_max)
        xknots_iteration.append(xknots.copy())
        yknots_iteration.append(yknots.copy())

        if j == 0:
            xknots_total = xknots.copy()
            yknots_total = yknots.copy()

            # otherwise, append the new knots to the total knots to build a smooth spline across the entire rc space
        else:
            # In the case of the x knots, simply append the new x knots to the total x knots less the first knot
            # of the total x knots list (this is repeated exactly in the new x knots list as its last value)
            xknots_total = np.append(xknots, xknots_total[1:])

            # In the case of the y_knot values:
            # 1. Do the same as the x knots AND
            # 2. Shift the new y knots by the first y knot value in the total knots list to avoid a discontinuity
            # in the spline fit at the junction of the two staircase steps
            yknots_total = np.append(yknots + yknots_total[0], yknots_total[1:])

    outer_fpt = fpt(xknots_total, yknots_total, xknots_total[0], xknots_total[-1], False)[0]

    return xknots_total, yknots_total, xknots_iteration, yknots_iteration, outer_fpt


def fpt_write_wrap_stair(name, stair_rc_list):
    print('input json:', f"{name}.json")
    with open(f"{name}.json", 'r') as input_json:
        data = json.load(input_json)

    set_defaults(data,
                 defaults={
                     "WL_sbias": 4.0,
                     "fail_max": 5,
                     "req_dists": 25000,
                     "rc_target_min_percentile": 0.10,
                     "nsteps_max": 1000,
                     "nboot": 0

                 })

    t_bonds = data['transient_bonds']
    #p_bonds = data['permanent_bonds']
    #nl_bonds = data['nonlocal_bonds']
    #max_d = data['req_dists']
    nboot = data['nboot']
    rc_transient = t_bonds[-1][-1]
    rh = data['rh']

    # Append the final step rc value and sort in descending order
    stair_rc_list.append(rc_transient)
    stair_rc_list.sort(reverse=1)

    t_ind, xb_on, yb_on, inner_fpt, t_on = inner_stair_fpt(name)
    x_total, y_total, x_iteration, y_iteration, outer_fpt = fpt_outer_stair(name, stair_rc_list)

    if plot_data or write_data:
        plot_outer_integrand(x_total, y_total, name)

    nknots = len(x_total)
    print(f"Completed fpts for {name} for outer nknots = {nknots - 1} and outer_fpt = {outer_fpt}.")

    output_i = [t_ind, inner_fpt, outer_fpt]
    if nboot:
        fpt_inner_std_val = fpt_std(t_on, rh, rc_transient, True, nboot, xb_on, yb_on)
        fpt_outer_std_val = fpt_outer_std(name, stair_rc_list, x_iteration, y_iteration, nboot)
        output_i += [fpt_inner_std_val, fpt_outer_std_val]
    else:
        output_i += [0.0, 0.0]
    output = []
    output.append(output_i)

    csv_name = f'{name}.csv'
    with open(csv_name, 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def inner_stair_fpt(name):
    with open(f"{name}.json", 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']

    t_bonds = data['transient_bonds']
    nl_bonds = data['nonlocal_bonds']
    max_d = data['req_dists']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    t_on = []

    dist_t_active = get_dists(name, t_ind, max_d)
    for j in range(len(dist_t_active)):
        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])

    t_on = np.array(t_on)
    xb, yb, fpt_on = fpt_per_bead_pair(t_on, rh, rc_transient, True, name)
    return t_ind, xb, yb, fpt_on, t_on


def knots_wrap(stair_output_name, rc_max, resample=False):
    with open(f"{stair_output_name}.json", 'r') as input_json:
        data = json.load(input_json)

    t_bonds = data['transient_bonds']
    nl_bonds = data['nonlocal_bonds']
    max_d = data['req_dists']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    dist_t_active = get_dists(stair_output_name, t_ind, max_d)

    t_on = []
    t_off = []

    for j in range(len(dist_t_active)):
        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])
        else:
            t_off.append(dist_t_active[j])

    t_off = np.array(t_off)

    if rc_max == -1.:
        rc_max = np.max(t_off)

    if (resample):
        boot_fpt_i = utils.resample(t_off)
        t_off = boot_fpt_i

    return find_knots(t_off, rc_transient, rc_max)


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


def integrate_spline(x_knot, y_knot):
    output = []
    f = lambda x: np.exp(-spline(x_knot, y_knot, x))

    for i in range(len(x_knot) - 1):
        output.append(integrate.quadrature(f, x_knot[i], x_knot[i + 1], maxiter=55, tol=1e-6, rtol=1e-6)[0])

    return sum(output)


def fun(x, y, x_i):
    #global best_args
    #global best_val

    integral = integrate_spline(x, y)

    f = np.sum(spline(x, y, x_i))
    count = len(x_i)

    if count > 0:
        f /= (count * 1.0)

    f += np.log(integral)

    #if f < best_val:
    #    best_val = f
    #    best_args = np.copy(y)

    # print('f = ', f, ' at y = ', y)
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
    if (i[0] == x_knot.size ):
        print('Error in dspline index')

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

    if BFGS:
        ub = np.full(len(x), 10.0)
        lb = np.full(len(x), -10.0)
        min_f = minimize(f, y0, method=nlopt.LD_LBFGS, jac=df, ub=ub, lb=lb)
    else:
        ub = np.full(len(x), 10.0)
        lb = np.full(len(x), -10.0)
        min_f = minimize(f, y0, method=nlopt.LN_COBYLA, jac=None, ub=ub, lb=lb)

    return min_f


def find_nearest(x_array, x_target):
    """
    Find the index of the closest value to x_target in x_array
    Parameters
    ----------
    x_array: np.array
        The array to search
    x_target: float

    Returns
    -------
    int
        The index of the closest value to x_target in x_array

    Examples
    --------
    Case 1: the target is in the array
    >>> x_array = np.array([ 1.8116996 ,  4.3631407 ,  4.72581521,  7.52534446,  9.78622537,
       10.18616866, 17.61249568, 20.07059104, 28.46182021, 29.39872047,
       35.45173343, 36.68139581, 46.64677056, 50.90701571, 56.85342466,
       62.56567633, 80.93815028, 92.66450544, 99.02348267, 99.96312029])
    >>> x_target = x_array[6]
    >>> find_nearest(x_array, x_target)
    5 # The x target 17.612495684496764 lies between 10.186168658934259 and 17.612495684496764

    Case 2: the target is not in the array
    >>> x_target = 20.0
    >>> find_nearest(x_array, x_target)
    6 # The x target 20.0 lies between 17.612495684496764 and 20.070591035714285
    """
    # find the index of the closest value to x_target in x less than x_target
    x_nearest_idx = np.searchsorted(x_array, x_target, side="left") - 1

    # print(f"The x target {x_target} lies between {x_array[x_nearest_idx]} and {x_array[x_nearest_idx + 1]}")
    return x_nearest_idx


def cdf_at_x(x_knot, y_knot, x, norm, cdf_base):
    """
    Function to calculate the cdf at a given x value
    Parameters
    ----------
    x_knot: np.array
    y_knot: np. Array
    x: float
    norm: float
    cdf_base: np.array

    Returns
    -------
    float
        The cdf at x
    """
    if x == x_knot[-1]:
        return 1.
    if x == x_knot[0]:
        return 0.

    x_nearest_idx = find_nearest(x_array=x_knot, x_target=x)

    pdf_func = lambda t: np.exp(-spline(x_knot, y_knot, t)) / norm

    pdf_integral = integrate.quadrature(pdf_func, x_knot[x_nearest_idx], x, maxiter=55, tol=1e-3, rtol=1e-3)[0]

    return cdf_base[x_nearest_idx] + pdf_integral


def fpt_integrand(x_knot, y_knot, x, state, norm):
    if type(state) != bool:
        raise TypeError('the variable state must be bool')

    free_energy = spline(x_knot, y_knot, x)

    pdf_func = lambda t: np.exp(-spline(x_knot, y_knot, t)) / norm
    pdf_val = np.exp(-free_energy) / norm

    # create a database of all the integrals over the knot intervals
    cdf_base = []
    # iterate through all the x knot values
    for i, el in enumerate(x_knot):
        # evaluate the integral from the first knot to the current knot of pdf_func
        res = integrate.quadrature(pdf_func, x_knot[0], el, tol=1e-3, rtol=1e-3,
                                   maxiter=90)[0]

        # print(f"The integral from {xknots[0]} to {el} is: {res}")
        # append this integral to the database
        cdf_base.append(res)

    # iterate through all the x values and find the cdf at each x value; accumulate the cdf values in a list "cdf"
    cdf = []
    for i in range(len(x)):
        cdf.append(cdf_at_x(x_knot, y_knot, x[i], norm, cdf_base))
    # convert the list to a numpy array
    cdf = np.array(cdf)

    if state is True:
        # fpt is inner
        integrand = cdf * cdf / pdf_val
    else:
        integrand = (1.0 - cdf) * (1.0 - cdf) / pdf_val

    return integrand


def plot_outer_integrand(x_knot, y_knot, name):
    norm = integrate_spline(x_knot, y_knot)
    x = np.linspace(x_knot[0], x_knot[-1], 1000)
    y_integrand = fpt_integrand(x_knot, y_knot, x, False, norm)
    pdf_val = np.exp(- spline(x_knot, y_knot, x)) / norm
    # bins = np.histogram_bin_edges(pdf_val, bins='auto', range=(x_knot[0], x_knot[-1]))
    # hist,bin_edges = np.histogram(pdf_val,density=True,bins='fd')

    # create a dictionary and pandas DataFrame
    my_dict = dict(x=x, y=pdf_val, z=y_integrand)
    data = pd.DataFrame(my_dict)

    if plot_data:
        fig, (ax1, ax2) = plt.subplots(1, 2)

        # sns.histplot(data=data,x='x',y='y',kde=True,ax=ax1)
        sns.lineplot(data=data, x='x', y='y', ax=ax1)

        sns.lineplot(data=data, x='x', y='z', ax=ax2)

        figName = f'{name}.png'

        plt.savefig(figName)
        plt.close()
    # plt.show()

    if write_data:
        datName = f'{name}.dat'
        file_object = open(datName, "w")
        for i in range(len(x)):
            print(x[i], pdf_val[i], y_integrand[i], file=file_object)
        file_object.close()


def fpt(x_knot, y_knot, xmin, xmax, state):
    f = lambda x: fpt_integrand(x_knot, y_knot, x, state, norm)
    norm = integrate_spline(x_knot, y_knot)
    fval = integrate.quadrature(f, xmin, xmax, tol=1e-3, rtol=1e-3, maxiter=100)
    return fval


def fpt_per_bead_pair_old(dist_vec, nknots, min_dist, max_dist, state):
    dist_vec.sort()

    x = np.linspace(min_dist, max_dist, nknots)
    # initial guess for y-coordinates of knots
    y = np.zeros(nknots)

    y = minimize_f(x, dist_vec, y)

    return fpt(x, y, min_dist, max_dist, state)[0]


def fpt_std(dist_vec, min_dist, max_dist, state, nboot, x, y):
    sample_size = int(len(dist_vec))

    boot_fpt = []
    for i in range(nboot):
        boot_fpt_i = utils.resample(dist_vec, random_state=i)
        fpt_i = fpt_boot(boot_fpt_i, min_dist, max_dist, state, x, y)
        if Verbose_Bootstrap:
            print('bootstrap round', i, ' has fpt = ', fpt_i)
        boot_fpt.append(fpt_i)

    return np.std(boot_fpt, ddof=1)


def chisquare_fit(x_knot, y_knot, dist_vec, norm):
    # assumes dist_vec is sorted and norm is calculated beforehand
    # Make discrete probability by integrating pdf over numK equal probability intervals
    # Perform chi^2 test_save on interval counts, assuming dof = the number of knot parameters
    #
    pdf_func = lambda t: np.exp(-spline(x_knot, y_knot, t)) / norm
    numK = 500
    numPoints = dist_vec.size
    if (numPoints < 5 * numK):
        numK = int(numPoints / 5)

    indices = np.linspace(0, numPoints, numK + 1, dtype=np.int64)
    Expected = np.zeros(numK)
    Observed = np.zeros(numK)
    for i in range(numK):
        first_index = indices[i]
        second_index = indices[i + 1]
        if i == numK - 1:
            second_index = numPoints - 1

        # distance values at evenly-spaced probability intervals of 1/numK
        x_i = dist_vec[first_index]
        x_f = dist_vec[second_index]

        Observed[i] = (second_index - first_index)
        Expected[i] = integrate.quadrature(pdf_func, x_i, x_f, tol=1e-3, rtol=1e-3, maxiter=95)[0] * (numPoints - 1)
        # print('range [', x_i, ',', x_f, ']. Expected[i] = ', Expected[i], ' percentiles = ', Observed[i])

    norm_p = np.sum(Observed)
    norm_e = np.sum(Expected)
    # must normalized expected results so that the counts agree to suitable precision for stats.chisquare
    Expected = Expected * norm_p / norm_e

    statistic, pval = stats.chisquare(f_obs=Observed, f_exp=Expected, ddof=x_knot.size)
    if Verbose_Convergence:
        print('      Chi-squared statistic: ', statistic, 'p-value:', pval, ' numK = ', numK, ' num knots = ',
              x_knot.size)
    # return probability that the statistic would arise by chance given the pdf
    return pval


def KuipersTest(x_knot, y_knot, dist_vec, norm):
    f = lambda t: np.exp(-spline(x_knot, y_knot, t))

    npoints = len(dist_vec)
    # print("In KStest with npoints = ", npoints)
    cdf = []
    cdf.append(0.)

    ecdfs = np.arange(npoints + 1, dtype=float) / npoints

    cdf_prev = 0.
    for i in range(npoints - 1):
        integral_val = integrate.quadrature(f, dist_vec[i], dist_vec[i + 1], maxiter=65, tol=1e-3, rtol=1e-3)[0] / norm
        cdf_i = cdf_prev + integral_val
        cdf.append(cdf_i)
        cdf_prev = cdf_i

    cdf.append(1.0)
    cdf = np.array(cdf)

    g = q = 0.0
    h1 = -np.inf
    h2 = -np.inf

    devIndex1 = -1
    devIndex2 = -1
    b = np.sqrt(npoints)
    for i in range(npoints):
        d1 = cdf[i] - ecdfs[i]

        if (d1 > h1):
            devIndex1 = i
            h1 = d1

        d2 = ecdfs[i + 1] - cdf[i]
        if (d2 > h2):
            devIndex2 = i
            h2 = d2

    h = h1 + h2
    d = b * h + 0.155 * h + 0.24 * h / b
    a = 2 * d * d

    for i in range(1, 201):
        term = (4.0 * a * i * i - 2.0) * np.exp(-a * i * i)
        q += term
        if np.abs(term) <= g:
            break
        else:
            g = np.abs(term) / 1000.

    maxDev = h
    devIndex = devIndex1
    if np.abs(d1) < np.abs(d2):
        devIndex = devIndex2

    devX = dist_vec[devIndex]

    # print(' h1 = ', h1, ' h2 = ', h2, ' for Kuipers test_save')
    if Verbose_Convergence:
        print('Kuipers test_save: q=', q, ' refining devX =', devX, ' num knots is ', x_knot.size, ' npoints = ',
              npoints)

    return q, maxDev, devX


def KStest(x_knot, y_knot, dist_vec, norm):
    f = lambda t: np.exp(-spline(x_knot, y_knot, t))

    npoints = len(dist_vec)
    # print("In KStest with npoints = ", npoints)
    cdf = []
    cdf.append(0.)

    ecdfs = np.arange(npoints + 1, dtype=float) / npoints
    # print(" ecdf is ", ecdfs)

    cdf_prev = 0.
    for i in range(npoints - 1):
        integral_val = integrate.quadrature(f, dist_vec[i], dist_vec[i + 1], maxiter=62, tol=1e-3, rtol=1e-3)[0] / norm
        cdf_i = cdf_prev + integral_val
        cdf.append(cdf_i)
        cdf_prev = cdf_i

    cdf.append(1.0)
    cdf = np.array(cdf)
    # print(' size of ecdf is ', ecdfs.size, ' size of cdf is ', cdf.size)
    # print('cdf is ', cdf)

    e = 0.
    devIndex = -1
    h = -1
    b = np.sqrt(npoints)
    c = 2.
    for i in range(npoints):
        d = np.fabs(cdf[i] - e)
        # print(' d is ', d, ' for i = ', i)
        if (d > h):
            devIndex = i
            h = d
        e = ecdfs[i + 1]
        d = np.fabs(e - cdf[i])

        if (d > h):
            devIndex = i
            h = d

    d = b * h + 0.12 * h + 0.11 * h / b
    q = kolmogorov(d)

    # print(' Calculate q value using h = ', h, ' Kn = ', Kn, ' devIndex = ', devIndex)
    maxDev = h
    devX = dist_vec[devIndex]
    if Verbose_Convergence:
        print('Adaptive KS: q =', q, ' devX = ', devX, ' num knots is ', x_knot.size)

    return q, maxDev, devX


def find_nearest_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def find_knots(dist_vec, min_dist, max_dist):
    # print('analyzing ', dist_vec.size, ' distances between ', min_dist, ' and ', max_dist)
    dist_vec.sort()
    if (dist_vec[0] < min_dist):
        print('Error: minimum is below what is allowed.')
    elif (dist_vec[-1] > max_dist):
        print('Error: maximum distance of ', dist_vec[-1], ' is above outer stair', max_dist)
    # delta_x = (max_dist - min_dist) / 20.
    q = 0.
    q_best = 0.
    nknots = 5
    x = np.linspace(min_dist, max_dist, nknots)
    y = np.zeros(nknots)

    while q < q_cut and x.size < 18:
        if x[0] != min_dist or x[-1] != max_dist:
            print('Error in initial values of xknot: x = ', x, ' for nknots = ', nknots, ' = ', x.size)

        y = minimize_f(x, dist_vec, y)
        norm = integrate_spline(x, y)

        if Adaptive:
            if Adaptive_KS:
                #  Use maximum deviation in KS test_save to position knots
                q, maxDev, devX = KStest(x, y, dist_vec, norm)
            else:
                # Use Kuipers test_save to refine knot positions:  supposedly more sensitive to tails of density
                q, maxDev, devX = KuipersTest(x, y, dist_vec, norm)

            if q > q_best:
                q_best = q
                x_best = copy.deepcopy(x)
                y_best = copy.deepcopy(y)

            i_nearest, x_nearest = find_nearest_value(x, devX)

            if x_nearest != devX:
                i_lower = find_nearest(x, devX)
                # print('devX = ', devX, ' is in interval [', x[i_lower], ',', x[i_lower + 1], ']')
                devX = x[i_lower] + 0.5 * (x[i_lower + 1] - x[i_lower])

            if q < q_cut:
                if Verbose_Convergence:
                    print('Will place knot at ', devX)

                x_new = copy.deepcopy(x)
                x_new = np.append(x_new, devX)
                x_new.sort()

                # check to see if devX is already one of the knots
                u, c = np.unique(x_new, return_counts=True)
                dup = u[c > 1]
                if len(dup) > 0:
                    if Verbose_Convergence:
                        print(' Adaptive process abandoned since suggested knot position ', devX,
                              ' duplicated in array ', x_new)

                    nknots = x.size
                    x = np.linspace(min_dist, max_dist, nknots)
                    y = np.zeros(nknots)
                else:
                    y_new = spline(x, y, x_new)
                    x = x_new
                    y = y_new
        else:
            q, maxDev, devX = KuipersTest(x, y, dist_vec, norm)
            if q > q_best:
                q_best = q
                x_best = copy.deepcopy(x)
                y_best = copy.deepcopy(y)

            if q < q_cut:
                x_new = np.linspace(min_dist, max_dist, nknots + 1)
                y_new = spline(x, y, x_new)
                x = x_new
                y = y_new
                nknots = x.size

    nknots = x_best.size
    q_chi = chisquare_fit(x_best, y_best, dist_vec, norm)
    if q_chi < 0.05:
        if Verbose_Convergence:
            print('Warning: chi-squared test_save value of q = ', q_chi, ' lies below generally accepted threshold.')
            print('    xmin = ', min_dist, ' xmax = ', max_dist)

    y_best -= y_best[-1]
    return nknots, x_best, y_best


def find_fixed_knots(dist_vec, x, y):
    # for fixed knot x values, find best y for data set dist_vec
    dist_vec.sort()
    y = minimize_f(x, dist_vec, y)
    y -= y[-1]
    return y


def fpt_boot(dist_vec, min_dist, max_dist, state, x, y):
    yi = find_fixed_knots(dist_vec, x, y)
    fpt_val = fpt(x, yi, min_dist, max_dist, state)
    return fpt_val[0]


def fpt_per_bead_pair(dist_vec, min_dist, max_dist, state, struc_id):
    dist_vec.sort()

    nknots, xb, yb = find_knots(dist_vec, min_dist, max_dist)
    fpt_val = fpt(xb, yb, min_dist, max_dist, state)[0]

    if state:
        print(f"{struc_id} inner fpt: Converged for nknots = {nknots}  with {len(dist_vec)} distances: fpt = {fpt_val}")
    else:
        print(f"{struc_id} outer fpt: Converged for nknots = {nknots}  with {len(dist_vec)} distances: fpt= {fpt_val}")

    if (plot_data or write_data) and state == False:
        plot_outer_integrand(xb, yb, struc_id)

    return xb, yb, fpt_val


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
