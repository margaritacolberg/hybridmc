import csv
import json
import os

import h5py
import nlopt
import numpy as np
from scipy import sparse, interpolate, integrate
from scipy.special import kolmogorov
from sklearn import utils
from sys import path

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

path.append('../')
from post_processing import matrix_element


# tol, rtol, maxiter = 1e-3, 1e-3, 100
plot_data = True
BFGS = True

best_val = np.inf
best_args = None

import numpy as np


def minimize(f, x0, method, jac, ub, lb):
    global best_args
    global best_val

    dim = len(x0)

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

    best_args = np.zeros(dim)
    best_val = np.inf

    try:
        opt_results = opt.optimize(x0)
    except:
        print("Failure: result code = ", opt.last_optimize_result())
        print('Trying non-derivative COBYLA routine.')
        print('opt_results are ', opt.last_optimum_value(), ' = ', best_val )
        print('best yvalues are ', best_args)
        opt2 = nlopt.opt(nlopt.LN_NELDERMEAD,dim)
        opt2.set_min_objective(obj_f)
        opt2.set_ftol_rel(1e-4)
        opt2.set_xtol_rel(1e-4)

        opt_results = opt2.optimize(best_args)

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


def fpt_write_wrap_default(name, nboot=0, boot_name=""):
    print('input json:', f"{name}.json")
    with open(f"{name}.json", 'r') as input_json:
        data = json.load(input_json)

    rh = data['rh']

    t_bonds = data['transient_bonds']
    p_bonds = data['permanent_bonds']
    nl_bonds = data['nonlocal_bonds']

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
    if nboot and name == boot_name:
        run_bootstrap = True
    else:
        run_bootstrap = False

    output = []

    # generate for non staircase simulations, the values of the distances within and outside the rc
    t_on = []
    t_off = []

    dist_t_active = get_dists(name, t_ind)
    for j in range(len(dist_t_active)):

        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])
        else:
            t_off.append(dist_t_active[j])

    t_on = np.array(t_on)
    t_off = np.array(t_off)

    # inner fpt for transient bond turned on
    fpt_on = fpt_per_bead_pair(t_on, rh, rc_transient, True, name)[0]
    fpt_off = fpt_per_bead_pair(t_off, rc_transient, np.max(t_off), False, name)[0]

    output_i = [t_ind, fpt_on, fpt_off]

    if run_bootstrap:
        fpt_on_var = fpt_var(t_on, rh, rc_transient, True, nboot, name)
        fpt_off_var = fpt_var(t_off, rc_transient, np.max(t_off), False, nboot, name)
        output_i += [fpt_on_var, fpt_off_var]

    output.append(output_i)

    with open(f"{name}.csv", 'w') as output_csv:
        writer = csv.writer(output_csv)
        writer.writerows(output)


def fpt_write_wrap_stair(name, stair_rc_list):
    print('input json:', f"{name}.json")
    with open(f"{name}.json", 'r') as input_json:
        data = json.load(input_json)

    t_bonds = data['transient_bonds']
    rc_transient = t_bonds[-1][-1]

    # Append the final step rc value and sort in descending order
    stair_rc_list.append(rc_transient)
    stair_rc_list.sort(reverse=1)

    # get knots for the whole rc space to make one smooth spline fit across the outer range encapsulating all the rcs
    xknots_total = []
    yknots_total = []
    for j in range(len(stair_rc_list)):

        # So long as the rc is not at the target we append the intermediate rc to the output name
        if stair_rc_list[j] != stair_rc_list[-1]:
            stair_output_name = f'{name}_{stair_rc_list[j]}'
        else:
            stair_output_name = name

        # if the rc is the outermost rc (step 1), then the maximum rc can be the larget distance in the data file
        # indicate this by using -1 as the rc_max
        if j == 0:
            rc_max = -1.

        # otherwise, the maximum rc is the rc of the previous step
        else:
            rc_max = stair_rc_list[j - 1]

        nknots, xknots, yknots = knots_wrap(stair_output_name, rc_max)

        # if this is the first step, then the knots are the total knots that subsequent knots are added
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

    t_ind, inner_fpt = inner_stair_fpt(name)

    outer_fpt = fpt(xknots_total, yknots_total, xknots_total[0], xknots_total[-1], False)[0]

    if plot_data:
        plot_outer_integrand(xknots_total, yknots_total, name)



    nknots = len(xknots_total)
    print(f"Completed fpts for {name} for outer nknots = {nknots - 1} and outer_fpt = {outer_fpt}.")

    output = []
    output.append([t_ind, inner_fpt, outer_fpt])
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

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    t_on = []

    dist_t_active = get_dists(name, t_ind)
    for j in range(len(dist_t_active)):
        if dist_t_active[j] < rc_transient:
            t_on.append(dist_t_active[j])

    t_on = np.array(t_on)
    fpt_on = fpt_per_bead_pair(t_on, rh, rc_transient, True, name)[0]
    return t_ind, fpt_on


def knots_wrap(stair_output_name, rc_max):
    with open(f"{stair_output_name}.json", 'r') as input_json:
        data = json.load(input_json)

    t_bonds = data['transient_bonds']
    nl_bonds = data['nonlocal_bonds']

    rc_transient = t_bonds[-1][-1]

    # get index of the transient bond in the list of nonlocal bonds
    t_ind = nl_bonds.index(t_bonds[0])

    dist_t_active = get_dists(stair_output_name, t_ind)

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
        output.append(integrate.quadrature(f, x_knot[i], x_knot[i + 1], tol=1e-6,
                                           rtol=1e-6)[0])

    return sum(output)


def fun(x, y, x_i):
    global best_args
    global best_val

    integral = integrate_spline(x, y)

    f = np.sum(spline(x, y, x_i))
    count = len(x_i)

    if count > 0:
        f /= (count * 1.0)

    f += np.log(integral)

    if f < best_val:
        best_val = f
        best_args = np.copy(y)

    #print('f = ', f, ' at y = ', y)
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


    if BFGS:
        min_f = minimize(f, y0, method=nlopt.LD_LBFGS, jac=df, ub=None, lb=None)
    else:
        ub = np.full(len(x),5.0)
        lb = np.full(len(x),-5.0)
        min_f = minimize(f, y0, method=nlopt.LN_COBYLA, jac=None, ub=ub,lb=lb)


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

    pdf_integral = integrate.quadrature(pdf_func, x_knot[x_nearest_idx], x, tol=1e-6, rtol=1e-6)[0]

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
        res = integrate.quadrature(pdf_func, x_knot[0], el, tol=1e-6, rtol=1e-6,
                                   maxiter=100)[0]

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



def plot_outer_integrand(x_knot,y_knot,name):
    norm = integrate_spline(x_knot, y_knot)
    x = np.linspace(x_knot[0],x_knot[-1],1000)
    y_integrand = fpt_integrand(x_knot, y_knot, x, False, norm)
    pdf_val = np.exp(- spline(x_knot, y_knot,x) ) / norm
    #bins = np.histogram_bin_edges(pdf_val, bins='auto', range=(x_knot[0], x_knot[-1]))
    #hist,bin_edges = np.histogram(pdf_val,density=True,bins='fd')

    # create a dictionary and pandas DataFrame
    my_dict = dict(x=x,y=pdf_val,z=y_integrand)
    data = pd.DataFrame (my_dict)

    fig, (ax1,ax2) = plt.subplots(1,2)

    #sns.histplot(data=data,x='x',y='y',kde=True,ax=ax1)
    sns.lineplot(data=data,x='x',y='y',ax=ax1)
    sns.lineplot(data=data,x='x',y='z',ax=ax2)

    figName = f'{name}.png'

    plt.savefig(figName)
    #plt.show()

    datName = f'{name}.dat'
    file_object = open(datName, "w")
    for i in range(len(x)):
        print(x[i], pdf_val[i],y_integrand[i], file=file_object)
    file_object.close()







def fpt(x_knot, y_knot, xmin, xmax, state):
    f = lambda x: fpt_integrand(x_knot, y_knot, x, state, norm)
    norm = integrate_spline(x_knot, y_knot)
    fval = integrate.quadrature(f, xmin, xmax, tol=1e-6, rtol=1e-6, maxiter=100)
    return fval


def fpt_per_bead_pair_old(dist_vec, nknots, min_dist, max_dist, state):
    dist_vec.sort()

    x = np.linspace(min_dist, max_dist, nknots)
    # initial guess for y-coordinates of knots
    y = np.zeros(nknots)

    y = minimize_f(x, dist_vec, y)

    return fpt(x, y, min_dist, max_dist, state)


def fpt_var(dist_vec, min_dist, max_dist, state, nboot, name):
    sample_size = int(len(dist_vec) / 2)

    boot_fpt = []
    for i in range(nboot):
        print('bootstrap round', i)
        boot_fpt_i = utils.resample(dist_vec, n_samples=sample_size,
                                    random_state=i)
        boot_fpt.append(fpt_per_bead_pair(boot_fpt_i, min_dist,
                                          max_dist, state, name)[0])

    return np.var(boot_fpt, ddof=1)


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
        integral_val = integrate.quadrature(f, dist_vec[i], dist_vec[i + 1], tol=1e-6, rtol=1e-6)[0] / norm
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
    Kn = b * h
    d = b * h + 0.12 * h + 0.11 * h / b
    q = kolmogorov(d)
    a = -2 * d * d
    g = 0.
    # print(' Calculate q value using h = ', h, ' Kn = ', Kn, ' devIndex = ', devIndex)
    maxDev = h
    # print('q is ', q, ' devIndex is ', devIndex, ' maxDev is ', maxDev, ' num knots is ', x_knot.size)
    devX = dist_vec[devIndex]

    return q, maxDev, devX


def find_knots(dist_vec, min_dist, max_dist):
    #print('analyzing ', dist_vec.size, ' distances between ', min_dist, ' and ', max_dist)
    dist_vec.sort()
    #print('dist_vec is ', dist_vec)

    q = 0.
    nknots = 4
    while q < 0.95 and nknots < 18:
        nknots = nknots + 1
        x = np.linspace(min_dist, max_dist, nknots)
        y = np.zeros(nknots)
        y = minimize_f(x, dist_vec, y)
        norm = integrate_spline(x, y)
        q, maxDev, devX = KStest(x, y, dist_vec, norm)
        #print('q is ', q, ' for nknots = ', nknots)

    # make last index zero
    y -= y[-1]

    return nknots, x, y


def fpt_per_bead_pair(dist_vec, min_dist, max_dist, state, struc_id):
    dist_vec.sort()

    q = 0
    nknots = 5
    while q < 0.95 and nknots < 18:
        x = np.linspace(min_dist, max_dist, nknots)
        # initial guess for y-coordinates of knots8
        y = np.zeros(nknots)

        y = minimize_f(x, dist_vec, y)

        norm = integrate_spline(x, y)
        # print('norm is ', norm)
        q, maxDev, devX = KStest(x, y, dist_vec, norm)
        nknots = nknots + 1

    if state:
        print(f"{struc_id} inner fpt: Converged for q = {q} for nknots = {nknots - 1}  with {len(dist_vec)} distances.")
    else:
        print(f"{struc_id} outer fpt: Converged for q = {q} for nknots = {nknots - 1}  with {len(dist_vec)} distances.")

    if plot_data and state == False:
        plot_outer_integrand(x,y,struc_id)

    return fpt(x, y, min_dist, max_dist, state)


def if_stair(ref_sim_id, files):
    """
    Function to check if the file_path has associated staircased steps simulations results. If yes, then the paths of these
    intermediate steps' csv files with their mfpt information are compiled into a list and returned.

    :param file_path: the final step mfpt simulation id
    :param files: the list of files in directory of interest
    :return: list[float] containing all the intermediate staircase rc values if they exist
    """

    # initialize output list merging all intermediate stair mfpt csv paths
    stair_rc_list = []

    # loop through json files in directory
    for file in files:
        if file.endswith('.json'):
            # obtain simulation tag for this file
            sim_id = file.rstrip('.json')
            # if the reference tag and this are the same then add the filepath to the output list
            if ref_sim_id.split('_') == sim_id.split('_')[:-1]:
                stair_rc_list.append(float(sim_id.split("_")[-1]))

    return stair_rc_list


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
