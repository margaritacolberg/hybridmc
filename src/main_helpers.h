//
// Created by vrajesh on 19/10/23.
//

#ifndef HYBRIDMC_MAIN_HELPERS_H
#define HYBRIDMC_MAIN_HELPERS_H

#include "crankshaft.h"
#include "entropy.h"
#include "hardspheres.h"
#include "json.hpp"
#include "snapshot.h"
#include "writer.h"
#include <H5Cpp.h>
#include <boost/program_options.hpp>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <numeric>

void initialize_pos(System &sys, Random &mt, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    std::optional<std::string> input_name,
                    std::optional<std::string> snapshot_name,
                    unsigned int t_bonds);

void initialize_system(System &sys, Random &mt, const Param &p, const Box &box,
                      UpdateConfig &update_config, Cells &cells,
                      EventQueue &event_queue);

void run_step(System &sys, const Param &p, const Box &box,
              UpdateConfig &update_config, CountBond &count_bond,
              Cells &cells, EventQueue &event_queue,
              unsigned int step, double del_t);

void run_trajectory_eq(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, CountBond &count_bond,
                       unsigned int iter,
                       DistWriter &dist_writer, std::vector<double> &dist);

void run_trajectory(System &sys, Random &mt, const Param &p, const Box &box,
                    std::vector<double> &dist, UpdateConfig &update_config,
                    UpdateConfigWriter &update_config_writer,
                    PosWriter &pos_writer, VelWriter &vel_writer,
                    ConfigWriter &config_writer, DistWriter &dist_writer,
                    std::set<Config> &store_config, ConfigInt &store_config_int,
                    CountBond &count_bond,
                    unsigned int iter);

Config run_trajectory_wl(System &sys, Random &mt, const Param &p,
                         const Box &box, UpdateConfig &update_config,
                         CountBond &count_bond,
                         unsigned int iter_wl,
                         bool record_dists = false,
                         std::vector<double>* dist = nullptr,
                         DistWriter* dist_writer = nullptr);

void wang_landau_process(System &sys, Random &mt, const Param &p, const Box &box,
                 UpdateConfig &update_config, CountBond &count_bond,
                 unsigned int nstates, std::vector<double> &s_bias,
                 DistWriter &dist_writer, std::vector<double> &dist);

void from_json(const nlohmann::json &json, Param &p);

#endif
