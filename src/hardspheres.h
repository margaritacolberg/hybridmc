// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYBRIDMC_HARDSPHERES_H
#define HYBRIDMC_HARDSPHERES_H

#include "box.h"
#include "cells.h"
#include "config.h"
#include "event.h"
#include "params.h"
#include "system.h"
#include "vec3.h"
#include <optional>
#include <random>
#include <vector>

#if NDEBUG
#define LOG_DEBUG(x)
#else
#define LOG_DEBUG(x)                                                           \
  std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl
#endif

extern double max_time;

using Random = std::mt19937;

void unit_sphere(Random &mt, double &x, double &y, double &z);

bool init_pos(std::vector<Vec3> &pos, const Box &box, Random &mt,
              const Param &p);
void draw_linear_chain(std::vector<Vec3> &pos, const Param &p);

void init_update_config(std::vector<Vec3> &pos, UpdateConfig &update_config,
                        const Box &box, const NonlocalBonds &transient_bonds);

void shift_vel_CM_to_0(std::vector<Vec3> &vel);

void init_vel(std::vector<Vec3> &vel, Random &mt, const double temp,
              const double m);

bool t_until_inner_coll(double xij, double yij, double zij, double vxij,
                        double vyij, double vzij, const double sigma2,
                        double &t);

void t_until_outer_coll(double xij, double yij, double zij, double vxij,
                        double vyij, double vzij, const double sigma2,
                        double &t);

bool v_after_coll(double xij, double yij, double zij, const double sigma2,
                  double &vxi, double &vyi, double &vzi, double &vxj,
                  double &vyj, double &vzj, const std::optional<double> &dU,
                  const double m);

double s_of_inner_event(const std::vector<double> &s_bias,
                        UpdateConfig update_config, const Config bond_mask);

double s_of_outer_event(const std::vector<double> &s_bias,
                        UpdateConfig update_config, const Config bond_mask);

void init_cells(const std::vector<Vec3> &pos, const Box &box, Cells &cells);

void delta_tpv(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
               const Box &box, const std::vector<double> &times, unsigned int i,
               unsigned int j, double &t, double &dx, double &dy, double &dz,
               double &dvx, double &dvy, double &dvz);

double get_rc2_inner(const std::tuple<Config, double> t_bond_mask_tuple,
                     const std::tuple<Config, double> p_bond_mask_tuple);

double get_rc2_outer(double rc2_inner, const Config t_bond_mask, std::optional<double> stair2,
                     UpdateConfig &update_config);

void if_coll(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
             double rh2, std::optional<double> stair2, const Box &box,
             const std::vector<uint64_t> &counter, EventQueue &event_queue,
             const std::vector<double> &times, unsigned int i, unsigned int j,
             const NonlocalBonds &transient_bonds,
             const NonlocalBonds &permanent_bonds, UpdateConfig &update_config,
             const unsigned int max_nbonds);

void iterate_coll(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
                  double rh2, std::optional<double> stair2, const Box &box,
                  const std::vector<uint64_t> &counter, EventQueue &event_queue,
                  const std::vector<double> &times, const Cells &cells,
                  unsigned int icell, unsigned int i,
                  const NonlocalBonds &transient_bonds,
                  const NonlocalBonds &permanent_bonds,
                  UpdateConfig &update_config, const unsigned int max_nbonds);

void walls_of_cell(unsigned int &zmin, unsigned int &zmax, unsigned int &ymin,
                   unsigned int &ymax, unsigned int &xmin, unsigned int &xmax,
                   const Cells &cells, unsigned int i);

void add_events_for_one_bead(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel, double rh2,
    std::optional<double> stair2, const Box &box, const std::vector<uint64_t> &counter,
    EventQueue &event_queue, const std::vector<double> &times,
    const Cells &cells, unsigned int i, const NonlocalBonds &transient_bonds,
    const NonlocalBonds &permanent_bonds, UpdateConfig &update_config,
    const unsigned int max_nbonds);

void add_events_for_all_beads(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
    unsigned int nbeads, double rh2, std::optional<double> stair2, const Box &box,
    const std::vector<uint64_t> &counter, EventQueue &event_queue,
    const std::vector<double> &times, const Cells &cells,
    const NonlocalBonds &transient_bonds, const NonlocalBonds &permanent_bonds,
    UpdateConfig &update_config, const unsigned int max_nbonds);

bool t_until_cell(double xi, double yi, double zi, const Box &box,
                  const Cells &cells, double vxi, double vyi, double vzi,
                  unsigned int ixo, unsigned int iyo, unsigned int izo,
                  double &dt, BeadCellEvent::Wall &wall, unsigned int &ixn,
                  unsigned int &iyn, unsigned int &izn);

void if_cell(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel, const Box &box,
             const std::vector<uint64_t> &counter, EventQueue &event_queue, const std::vector<double> &times,
             unsigned int i, const Cells &cells);

void init_cell_events(const std::vector<Vec3> &pos,
                      const std::vector<Vec3> &vel, unsigned int nbeads,
                      const Box &box, const std::vector<uint64_t> &counter,
                      EventQueue &event_queue, const std::vector<double> &times,
                      const Cells &cells);

void move_to_new_cell(Cells &cells, unsigned int i, unsigned int ixn,
                      unsigned int iyn, unsigned int izn);

void add_events_for_bead_after_crossing(
    const std::vector<Vec3> &pos, const std::vector<Vec3> &vel, double rh2,
    std::optional<double> stair2, const Box &box, const std::vector<uint64_t> &counter,
    EventQueue &event_queue, const std::vector<double> &times,
    const Cells &cells, unsigned int i, BeadCellEvent::Wall wall,
    const NonlocalBonds &transient_bonds, const NonlocalBonds &permanent_bonds,
    UpdateConfig &update_config, const unsigned int max_nbonds);

void if_nearest_bond(const std::vector<Vec3> &pos, const std::vector<Vec3> &vel,
                     const Box &box, const std::vector<uint64_t> &counter,
                     EventQueue &event_queue, const std::vector<double> &times,
                     unsigned int i, unsigned int j, double near_min2,
                     double near_max2);

void if_nnearest_bond(const std::vector<Vec3> &pos,
                      const std::vector<Vec3> &vel, const Box &box,
                      const std::vector<uint64_t> &counter,
                      EventQueue &event_queue, const std::vector<double> &times,
                      unsigned int i, unsigned int j, double nnear_min2,
                      double nnear_max2);

void init_nearest_bond_events(const std::vector<Vec3> &pos,
                              const std::vector<Vec3> &vel, unsigned int nbeads,
                              const Box &box,
                              const std::vector<uint64_t> &counter,
                              EventQueue &event_queue,
                              const std::vector<double> &times,
                              double near_min2, double near_max2);

void init_nnearest_bond_events(const std::vector<Vec3> &pos,
                               const std::vector<Vec3> &vel,
                               unsigned int nbeads, const Box &box,
                               const std::vector<uint64_t> &counter,
                               EventQueue &event_queue,
                               const std::vector<double> &times,
                               double nnear_min2, double nnear_max2);

bool check_overlap(const unsigned int i, const std::vector<Vec3> &pos,
                   const std::vector<Vec3> &vel,
                   const std::vector<double> &times, const double sigma2,
                   const Box &box);

bool check_bond(const double min_dist2, const double max_dist2, double dx,
                double dy, double dz);

bool check_local_dist(const std::vector<Vec3> &pos_trial, const Box &box,
                      const double near_min2, const double near_max2,
                      const double nnear_min2, const double nnear_max2);

void update_pos(Vec3 &pos, Vec3 &vel, double &time, double update_time);

bool final_vel(const std::vector<Vec3> &pos, std::vector<Vec3> &vel,
               const double sigma2, unsigned int i, unsigned int j,
               const Box &box, const std::optional<double> &dS, const double m);

double total_k_E(const std::vector<Vec3> &vel, const double m);

double compute_hamiltonian(const std::vector<Vec3> &vel,
                           const std::vector<double> &s_bias,
                           const Config &config, const double m);

void dist_between_nonlocal_beads(const std::vector<Vec3> &pos, const Box &box,
                                 const NonlocalBonds &nonlocal_bonds,
                                 std::vector<double> &dist);

bool process_event(const MinNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MaxNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MinNNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MaxNNearestEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MinNonlocalInnerEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MaxNonlocalInnerEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const MaxNonlocalOuterEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

bool process_event(const BeadCellEvent &ev, System &sys, const Param &p,
                   const Box &box, EventQueue &event_queue, Cells &cells,
                   UpdateConfig &update_config, CountBond &count_bond);

#endif
