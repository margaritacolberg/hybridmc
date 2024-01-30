import copy

import numpy as np
import nlopt
from scipy import sparse, interpolate, integrate
from scipy import stats
from scipy.special import kolmogorov
import bisect
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

plot_data = True
write_data = True
BFGS = True
chisquared_gof = True
q_cut = 0.2
Adaptive_KS = True

best_val = np.inf
best_args = None



numPoints = 2001
numK = 100
mu, sigma = 10.0, 1 # mean and standard deviation
gamma = 1.5

#  Make a complicated probability density that is difficult to fit
dist1 = np.random.normal(mu, sigma, numPoints)
dist2 = mu-4*sigma + gamma*np.random.standard_cauchy(numPoints)
dist_vec = np.append(dist1, dist2)
dist_vec = np.sort(dist_vec)

dist_vec = dist_vec[(dist_vec > -5) & (dist_vec < 15)]
x_min = np.min(dist_vec)
x_max = np.max(dist_vec)




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
        ub = np.full(len(x),10.0)
        lb = np.full(len(x),-10.0)
        min_f = minimize(f, y0, method=nlopt.LD_LBFGS, jac=df, ub=None, lb=None)
    else:
        ub = np.full(len(x),10.0)
        lb = np.full(len(x),-10.0)
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
    h = -np.inf
    b = np.sqrt(npoints)

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

    maxDev = h
    devX = dist_vec[devIndex]

    print('KS test: q=', q, ' refining devX =', devX, ' num knots is ', x_knot.size)
    return q, maxDev, devX


def KuipersTest(x_knot, y_knot, dist_vec, norm):
    f = lambda t: np.exp(-spline(x_knot, y_knot, t))

    npoints = len(dist_vec)
    # print("In KStest with npoints = ", npoints)
    cdf = []
    cdf.append(0.)

    ecdfs = np.arange(npoints + 1, dtype=float) / npoints

    cdf_prev = 0.
    for i in range(npoints - 1):
        integral_val = integrate.quadrature(f, dist_vec[i], dist_vec[i + 1], tol=1e-6, rtol=1e-6)[0] / norm
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
            h1= d1

        d2 = ecdfs[i+1] - cdf[i]
        if (d2 > h2):
            devIndex2 = i
            h2 = d2

    h = h1 + h2
    d = b * h + 0.155 * h + 0.24 * h / b
    a = 2 * d * d

    for i in range(1,201):
        term = (4.0 * a * i * i - 2.0) * np.exp(-a * i * i)
        q += term
        if np.abs(term) <= g:
            break
        else:
            g = np.abs(term)/1000.

    maxDev = h
    devIndex = devIndex1
    if np.abs(d1) < np.abs(d2):
        devIndex = devIndex2

    devX = dist_vec[devIndex]

    #print(' h1 = ', h1, ' h2 = ', h2, ' for Kuipers test')
    print('Kuipers test: q=', q, ' refining devX =', devX, ' num knots is ', x_knot.size)

    return q, maxDev, devX


def chisquare_fit(x_knot, y_knot, dist_vec, norm):
    # assumes dist_vec is sorted and norm is calculated beforehand
    # Make discrete probability by integrating pdf over numK equal probability intervals
    # Perform chi^2 test on interval counts, assuming dof = the number of knot parameters
    #
    pdf_func = lambda t: np.exp(-spline(x_knot, y_knot, t)) / norm
    numK = 500
    numPoints = dist_vec.size
    if (numPoints < 5*numK):
        numK = int(numPoints/5)

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
        Expected[i] = integrate.quadrature(pdf_func, x_i, x_f, tol=1e-6, rtol=1e-6,maxiter=100)[0] * (numPoints-1)
        #print('range [', x_i, ',', x_f, ']. Expected[i] = ', Expected[i], ' percentiles = ', Observed[i])

    norm_p = np.sum(Observed)
    norm_e = np.sum(Expected)
    # must normalized expected results so that the counts agree to suitable precision for stats.chisquare
    Expected = Expected * norm_p / norm_e


    # the dof in the chi-squared statistic is nknots -1 since normalization imposes one constraint
    statistic, pval = stats.chisquare(f_obs=Observed, f_exp=Expected, ddof=x_knot.size-1)
    print('Chi-squared statistic:', statistic, 'p-value:', pval, ' numK = ',numK, ' num knots = ', x_knot.size)
    # return probability that the statistic would arise by chance given the pdf
    return pval


def find_nearest_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def find_knots(dist_vec, min_dist, max_dist):
    #print('analyzing ', dist_vec.size, ' distances between ', min_dist, ' and ', max_dist)
    dist_vec.sort()

    delta_x = (max_dist - min_dist)/20.
    q = 0.
    q_best = 0.
    nknots = 5
    x = np.linspace(min_dist, max_dist,nknots)
    y = np.zeros(nknots)


    while q < q_cut and x.size < 18:

        y = minimize_f(x, dist_vec, y)
        norm = integrate_spline(x, y)

        if Adaptive_KS:
            #  Use maximum deviation in KS test to position knots
            q, maxDev, devX = KStest(x,y,dist_vec,norm)
        else:
            # Use Kuipers test to refine knot positions:  supposedly more sensitive to tails of density
            q, maxDev, devX = KuipersTest(x, y, dist_vec, norm)

        if q > q_best:
            q_best = q
            x_best = copy.deepcopy(x)
            y_best = copy.deepcopy(y)

        i_nearest, x_nearest = find_nearest_value(x,devX)

        if x_nearest != devX:
            i_lower = find_nearest(x,devX)
            print('devX = ', devX, ' is in interval [',x[i_lower],',',x[i_lower+1],']')

            interval = x[i_lower+1] - x[i_lower]
            if abs(devX - x_nearest) < 0.01*interval:
                print(' devX is too near a knot so placing at mid-interval')
                devX = x[i_lower] + 0.5*(x[i_lower+1] - x[i_lower])

            print('Will place knot at ', devX)

        if q < q_cut:

            x_new = copy.deepcopy(x)
            x_new = np.append(x_new,devX)
            x_new.sort()

            # check to see if devX is already one of the knots
            u, c = np.unique(x_new, return_counts=True)
            dup = u[c > 1]
            if len(dup)>0:
                print(' Adaptive process abandoned since suggested knot position ', devX, ' duplicated in array ', x_new)
                nknots = x.size
                x = np.linspace(min_dist, max_dist, nknots)
                y = np.zeros(nknots)
            else:
                y_new = spline(x,y,x_new)
                x = x_new
                y= y_new

    nknots = x_best.size
    q_chi = chisquare_fit(x_best, y_best, dist_vec, norm)
    return nknots, x_best, y_best

def plot_data(x_knot,y_knot,dist_vec):
    name = 'fit'
    norm = integrate_spline(x_knot, y_knot)
    x = np.linspace(x_knot[0],x_knot[-1],1000)
    pdf_val = np.exp(- spline(x_knot, y_knot,x) ) / norm

    # create a dictionary and pandas DataFrame
    my_dict = dict(x=x,y=pdf_val)
    data = pd.DataFrame (my_dict)

    if plot_data:
        fig, ax1 = plt.subplots(1,1)
        sns.kdeplot(data=dist_vec)
        #sns.histplot(data=data,x='x',y='dist_vec',kde=True,ax=ax1)
        nbins = int(dist_vec.size/50)
        ax1.hist(dist_vec,density=True,bins=nbins)
        sns.lineplot(data=data,x='x',y='y',ax=ax1)

        yvals = np.exp(-spline(x_knot,y_knot,x_knot))/norm
        plt.scatter(x_knot,yvals)

        figName = f'{name}.png'

        plt.savefig(figName)
        plt.close()
    #plt.show()

    if write_data:
        datName = f'{name}.dat'
        file_object = open(datName, "w")
        for i in range(len(x)):
            print(x[i], pdf_val[i],file=file_object)
        file_object.close()


nknots,xk,yk = find_knots(dist_vec, x_min, x_max)
print(' Optimal ', nknots, ' knots: x = ', xk, ' y = ', yk)
plot_data(xk,yk,dist_vec)