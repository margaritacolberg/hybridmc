// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYBRIDMC_CRANKSHAFT_H
#define HYBRIDMC_CRANKSHAFT_H

#include "hardspheres.h"

void swapMC(System &sys, Random &mt, const Box& box, const Param &p);

void rodrigues_rotation(const std::vector<Vec3> &pos, const double theta,
                        std::vector<Vec3> &pos_trial, const unsigned int ind,
                        const Box &box);

bool check_bond_if_crankshaft(const double min_dist2, const double max_dist2,
                              double dx, double dy, double dz);

bool check_local_dist_if_crankshaft(const std::vector<Vec3> &pos_trial,
                                    const Box &box, const double near_min2,
                                    const double near_max2,
                                    const double nnear_min2,
                                    const double nnear_max2);

bool check_nonlocal_dist(const std::vector<Vec3> &pos_trial, const Box &box,
                         const double rh2, const std::optional<double> stair2,
                         const NonlocalBonds &transient_bonds, const NonlocalBonds &permanent_bonds);

UpdateConfig config_int(const std::vector<Vec3> &pos_trial, const Box &box,
                        const NonlocalBonds &transient_bonds);

bool accept_move(const std::vector<double> &s_bias, UpdateConfig &orig_config,
                 UpdateConfig &trial_config, Random &mt);

void crankshaft(std::vector<Vec3> &pos,
                UpdateConfig &update_config,
                const Box &box,
                const double near_min2, const double near_max2,
                const double nnear_min2, const double nnear_max2,
                const double rh2,
                const std::optional<double> stair2,
                const NonlocalBonds &transient_bonds,
                const NonlocalBonds &permanent_bonds, Random &mt,
                const std::vector<double> &s_bias);

#endif
