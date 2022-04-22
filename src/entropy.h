// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYBRIDMC_ENTROPY_H
#define HYBRIDMC_ENTROPY_H

#include "hardspheres.h"

void init_s(std::vector<double> &s_bias, const unsigned int nbonds);

double compute_norm(const std::vector<uint64_t> &config_count,
                    const unsigned int nstates);

void compute_entropy(std::vector<uint64_t> &config_count,
                     std::vector<double> &s_bias, const unsigned int nstates,
                     const double pos_scale, const double neg_scale);

double compute_g_val(std::vector<uint64_t> &config_count,
                     const unsigned int nstates);

bool g_test(std::vector<uint64_t> &config_count, const unsigned int nstates,
            const double sig_level);

#endif
