// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYBRIDMC_SNAPSHOT_H
#define HYBRIDMC_SNAPSHOT_H

#include "hardspheres.h"
#include <string>

void write_snapshot(const std::string snapshot_name,
                    const std::vector<Vec3> &pos,
                    const std::vector<double> &s_bias, Random &mt,
                    UpdateConfig &update_config);

void write_ensemble(H5::H5Location &file, const System &sys);



void read_snapshot(const std::string snapshot_name, std::vector<Vec3> &pos,
                   std::vector<double> &s_bias, Random &mt,
                   UpdateConfig &update_config);

void read_input(const std::string input_name, std::vector<Vec3> &pos);

#endif
