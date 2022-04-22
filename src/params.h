// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYBRIDMC_PARAMS_H
#define HYBRIDMC_PARAMS_H

#include "config.h"
#include <optional>
#include <vector>

struct Param {
  // mass of each bead
  double m;
  // diameter of each bead
  double sigma;
  // squared diameter
  double sigma2;
  // distance constraints for nearest neighbors
  double near_min;
  double near_max;
  // smallest bond length squared between nearest neighbors
  double near_min2;
  // largest bond length squared between nearest neighbors
  double near_max2;
  // distance constrains for next-nearest neighbors
  double nnear_min;
  double nnear_max;
  // smallest bond length squared between next-nearest neighbors
  double nnear_min2;
  // largest bond length squared between next-nearest neighbors
  double nnear_max2;
  // distance constraints for nonlocal beads
  double rh;
  double rc;
  // smallest bond length squared between nonlocal beads
  double rh2;
  // largest bond length squared between nonlocal beads
  double rc2;
  // edge of previous stair in staircase potential of transient bond
  std::optional<double> stair;
  std::optional<double> stair2;
  // largest bond length between nonlocal beads in a permanent bond, if the
  // transient bond has a staircase potential
  std::optional<double> p_rc;
  std::optional<double> p_rc2;

  // vector of indices of beads which form bonds
  NonlocalBonds nonlocal_bonds;
  // vector of indices of beads which form bonds that can be broken
  NonlocalBonds transient_bonds;
  // vector of indices of beads which form bonds that cannot be broken
  NonlocalBonds permanent_bonds;
  // vector of indices of beads whose potential is a staircase
  std::optional<NonlocalBonds> stair_bonds;

  // number of tries to place beads in box
  uint64_t tries;
  // number of beads
  unsigned int nbeads;
  // length of box
  double length;
  // number of cells in each dimension of the box
  unsigned int ncell;
  // number of time intervals
  unsigned int nsteps;
  // step at which to write output to file
  unsigned int write_step;
  // increment time
  double del_t;
  // seeds for random number generator
  std::vector<unsigned int> seeds;

  // temperature of the system
  double temp;
  // number of Monte Carlo crankshaft moves
  unsigned int mc_moves;
  // number of times MC crankshaft moves are carried out before configuration
  // saved
  unsigned int mc_write;
  // total number of iterations of MD trajectories and MC moves
  unsigned int total_iter;
  // entropy scaling factors
  double pos_scale;
  double neg_scale;
  // probability that represents region of a distribution where results are
  // statistically significant
  double sig_level;
  // max number of G-test executions
  unsigned int max_g_test_count;
  // max number of nonlocal bonds that can form
  unsigned int max_nbonds;
};

#endif
