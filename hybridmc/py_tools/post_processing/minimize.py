#!/usr/bin/env python3
# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause
#
# minimize.py determines the minimum of a function using nonlinear optimization
#
# NLopt library: https://nlopt.readthedocs.io/en/latest/NLopt_Python_Reference/

import nlopt


def minimize(f, x0, method, jac, ub, lb):
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

    return opt.optimize(x0)


def make_nlopt_f(f, jac):
    def nlopt_f(x, grad):

        if grad.size > 0:
            grad[:] = jac(x)

        return f(x)
    return nlopt_f
