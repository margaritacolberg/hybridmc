// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// hybridmc_test.cc contains unit tests for functions contained in
// hardspheres.cc, entropy.cc, snapshot.cc, etc.

#include "cells.h"
#include "config.h"
#include "crankshaft.h"
#include "entropy.h"
#include "hardspheres.h"
#include "snapshot.h"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardspheres
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(inner_collision_time) {
  const double sigma2 = 1.0 * 1.0;

  // particles moving away from each other along the x-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(2.0, 0.0, 0.0, 10.0, 0.0, 0.0, sigma2, del_t));
  }
  // particles moving away from each other along the y-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(0.0, 3.0, 0.0, 0.0, 5.0, 0.0, sigma2, del_t));
  }
  // particles moving away from each other along the z-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(0.0, 0.0, 5.0, 0.0, 0.0, 8.0, sigma2, del_t));
  }
  // particles moving toward each other along the x-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        t_until_inner_coll(2.0, 0.0, 0.0, -10.0, 0.0, 0.0, sigma2, del_t));
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.1, 1e-10);
  }
  // particles moving toward each other along the y-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        t_until_inner_coll(0.0, 3.0, 0.0, 0.0, -5.0, 0.0, sigma2, del_t));
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.4, 1e-10);
  }
  // particles moving toward each other along the z-dimension
  {
    double del_t = 0;
    BOOST_CHECK(
        t_until_inner_coll(0.0, 0.0, 5.0, 0.0, 0.0, -8.0, sigma2, del_t));
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.5, 1e-10);
  }
  // particles moving toward each other along x but are offset in y by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(2.0, 3.0, 0.0, -10.0, 0.0, 0.0, sigma2, del_t));
  }
  // particles moving toward each other along x but are offset in z by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(2.0, 0.0, 5.0, -10.0, 0.0, 0.0, sigma2, del_t));
  }
  // particles moving toward each other along y but are offset in x by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(2.0, 3.0, 0.0, 0.0, -5.0, 0.0, sigma2, del_t));
  }
  // particles moving toward each other along y but are offset in z by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(0.0, 3.0, 5.0, 0.0, -5.0, 0.0, sigma2, del_t));
  }
  // particles moving toward each other along z but are offset in x by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(2.0, 0.0, 5.0, 0.0, 0.0, -8.0, sigma2, del_t));
  }
  // particles moving toward each other along z but are offset in y by more
  // than sigma
  {
    double del_t = 0;
    BOOST_CHECK(
        !t_until_inner_coll(0.0, 3.0, 5.0, 0.0, 0.0, -8.0, sigma2, del_t));
  }
}

BOOST_AUTO_TEST_CASE(outer_collision_time) {
  const double sigma2 = 1.17 * 1.17;

  // particles moving toward each other along the x-dimension; output should
  // equal (sigma + dx) / dvx
  {
    double del_t = 0;
    t_until_outer_coll(0.2, 0.0, 0.0, -10.0, 0.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.137, 1e-10);
  }
  // particles moving toward each other along the y-dimension
  {
    double del_t = 0;
    t_until_outer_coll(0.0, 0.3, 0.0, 0.0, -5.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.294, 1e-10);
  }
  // particles moving toward each other along the z-dimension
  {
    double del_t = 0;
    t_until_outer_coll(0.0, 0.0, 0.5, 0.0, 0.0, -8.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.20875, 1e-10);
  }
  // particles moving away from each other along the x-dimension
  {
    double del_t = 0;
    t_until_outer_coll(0.2, 0.0, 0.0, 10.0, 0.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.097, 1e-10);
  }
  // particles moving away from each other along the y-dimension
  {
    double del_t = 0;
    t_until_outer_coll(0.0, 0.3, 0.0, 0.0, 5.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.174, 1e-10);
  }
  // particles moving away from each other along the z-dimension
  {
    double del_t = 0;
    t_until_outer_coll(0.0, 0.0, 0.5, 0.0, 0.0, 8.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.08375, 1e-10);
  }
  // particles moving away from each other along x but are offset in y by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(0.2, 3.0, 0.0, 10.0, 0.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
  // particles moving away from each other along x but are offset in z by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(0.2, 0.0, 5.0, 10.0, 0.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
  // particles moving away from each other along y but are offset in x by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(2.0, 0.3, 0.0, 0.0, 5.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
  // particles moving away from each other along y but are offset in z by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(0.0, 0.3, 5.0, 0.0, 5.0, 0.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
  // particles moving away from each other along z but are offset in x by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(2.0, 0.0, 0.5, 0.0, 0.0, 8.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
  // particles moving away from each other along z but are offset in y by
  // more than sigma
  {
    double del_t = 0;
    t_until_inner_coll(0.0, 3.0, 0.5, 0.0, 0.0, 8.0, sigma2, del_t);
    BOOST_CHECK_CLOSE_FRACTION(del_t, 0.0, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(init_velocities) {
  const unsigned int nbeads = 10;
  std::vector<Vec3> vel(nbeads);
  const double temp = 1.5, m = 3.0;

  const unsigned int seed = 40;
  Random mt(seed);

  init_vel(vel, mt, temp, m);

  double vel_x_tot = 0.0, vel_y_tot = 0.0, vel_z_tot = 0.0;
  for (unsigned int i = 0; i < nbeads; i++) {
    vel_x_tot += vel[i].x;
    vel_y_tot += vel[i].y;
    vel_z_tot += vel[i].z;
  }

  // center of mass velocity
  double vel_x_CM = vel_x_tot / double(nbeads);
  double vel_y_CM = vel_y_tot / double(nbeads);
  double vel_z_CM = vel_z_tot / double(nbeads);

  BOOST_CHECK_SMALL(vel_x_CM, 1e-15);
  BOOST_CHECK_SMALL(vel_y_CM, 1e-15);
  BOOST_CHECK_SMALL(vel_z_CM, 1e-15);
}

BOOST_AUTO_TEST_CASE(collision_velocities) {
  // reflection at sigma, dS = {}
  {
    double vxi = -10.0;
    double vyi = 0.0;
    double vzi = 0.0;
    double vxj = 10.0;
    double vyj = 0.0;
    double vzj = 0.0;
    v_after_coll(-1.0, 0.0, 0.0, 1.0, vxi, vyi, vzi, vxj, vyj, vzj, {}, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(vxi, 10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vxj, -10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyj, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzj, 0.0, 1e-10);
  }

  // reflection at sigma, dS = {}
  {
    double vxi = 10.0;
    double vyi = 0.0;
    double vzi = 0.0;
    double vxj = -10.0;
    double vyj = 0.0;
    double vzj = 0.0;
    v_after_coll(1.0, 0.0, 0.0, 1.0, vxi, vyi, vzi, vxj, vyj, vzj, {}, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(vxi, -10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vxj, 10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyj, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzj, 0.0, 1e-10);
  }

  // transmission at rc with particles entering well
  {
    double vxi = -10.0;
    double vyi = 0.0;
    double vzi = 0.0;
    double vxj = 10.0;
    double vyj = 0.0;
    double vzj = 0.0;
    v_after_coll(2.25, 0.0, 0.0, 2.25, vxi, vyi, vzi, vxj, vyj, vzj, 0.0, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(vxi, -10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vxj, 10.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyj, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzj, 0.0, 1e-10);
  }

  // transmission at rc with particles entering well
  {
    double vxi = 10.0;
    double vyi = 0.0;
    double vzi = 0.0;
    double vxj = -10.0;
    double vyj = 0.0;
    double vzj = 0.0;
    v_after_coll(2.25, 0.0, 0.0, 2.25, vxi, vyi, vzi, vxj, vyj, vzj, -5.0, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(vxi, 10.2486263321547, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vxj, -10.2486263321547, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyj, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzj, 0.0, 1e-10);
  }

  // transmission at rc with particles exiting well
  {
    double vxi = -10.0;
    double vyi = 0.0;
    double vzi = 0.0;
    double vxj = 10.0;
    double vyj = 0.0;
    double vzj = 0.0;
    v_after_coll(2.25, 0.0, 0.0, 2.25, vxi, vyi, vzi, vxj, vyj, vzj, 5.0, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(vxi, -9.74859546128699, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzi, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vxj, 9.74859546128699, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vyj, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(vzj, 0.0, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(boundary_conditions) {
  const Box box{10};
  {
    double x = 12.0;
    double y = 18.0;
    double z = 15.0;
    box.minpos(x, y, z);
    BOOST_CHECK_CLOSE_FRACTION(x, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(y, 8.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(z, 5.0, 1e-10);
  }

  {
    double x = -198.0;
    double y = -356.0;
    double z = -221.0;
    box.minpos(x, y, z);
    BOOST_CHECK_CLOSE_FRACTION(x, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(y, 4.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(z, 9.0, 1e-10);
  }

  {
    double x = 3.0;
    double y = 7.0;
    double z = 6.0;
    box.minpos(x, y, z);
    BOOST_CHECK_CLOSE_FRACTION(x, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(y, 7.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(z, 6.0, 1e-10);
  }

  {
    double dx = 12.0;
    double dy = 18.0;
    double dz = 15.0;
    box.mindist(dx, dy, dz);
    BOOST_CHECK_CLOSE_FRACTION(dx, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dy, -2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dz, -5.0, 1e-10);
  }

  {
    double dx = -198.0;
    double dy = -356.0;
    double dz = -221.0;
    box.mindist(dx, dy, dz);
    BOOST_CHECK_CLOSE_FRACTION(dx, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dy, 4.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dz, -1.0, 1e-10);
  }

  {
    double dx = 3.0;
    double dy = 7.0;
    double dz = 6.0;
    box.mindist(dx, dy, dz);
    BOOST_CHECK_CLOSE_FRACTION(dx, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dy, -3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(dz, -4.0, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(cell_position) {
  const unsigned int ncell = 3;
  const double lcell = 1.0;

  const std::vector<Vec3> pos = {
      {1.8, 0.5, 2.3}, {0.0, 0.5, 0.0}, {2.0, 2.3, 1.4}, {0.5, 0.5, 2.7},
      {2.0, 2.0, 2.0}, {0.5, 1.8, 0.2}, {2.2, 2.6, 1.9}, {1.5, 2.0, 3.0},
      {1.0, 2.7, 0.6}, {2.5, 1.5, 2.5}};

  Cells cells(ncell, lcell);
  init_cells(pos, ncell, cells);

  BOOST_REQUIRE_EQUAL(cells.cell_x.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_y.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_z.size(), 10);
  BOOST_CHECK_EQUAL(cells.cell_x[0], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[0], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[0], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[2], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[2], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[2], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[3], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[3], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[3], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[5], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[6], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[6], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[6], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[7], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[7], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[7], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[8], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[8], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[8], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[9], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[9], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[9], 2);

  BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 1);
  // first index is 1D index of cell, second index is index of particle in
  // list of particles contained within cell, third value (in this case, 1)
  // is the index of particle within the pos vector
  BOOST_CHECK_EQUAL(cells.cells[0][0], 1);

  BOOST_REQUIRE_EQUAL(cells.cells[1].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[2].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[3].size(), 1);
  BOOST_CHECK_EQUAL(cells.cells[3][0], 5);

  BOOST_REQUIRE_EQUAL(cells.cells[4].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[5].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[6].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[7].size(), 2);
  BOOST_CHECK_EQUAL(cells.cells[7][0], 7);
  BOOST_CHECK_EQUAL(cells.cells[7][1], 8);

  BOOST_REQUIRE_EQUAL(cells.cells[8].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[9].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[10].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[11].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[12].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[13].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[14].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[15].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[16].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[17].size(), 2);
  BOOST_CHECK_EQUAL(cells.cells[17][0], 2);
  BOOST_CHECK_EQUAL(cells.cells[17][1], 6);

  BOOST_REQUIRE_EQUAL(cells.cells[18].size(), 1);
  BOOST_CHECK_EQUAL(cells.cells[18][0], 3);

  BOOST_REQUIRE_EQUAL(cells.cells[19].size(), 1);
  BOOST_CHECK_EQUAL(cells.cells[19][0], 0);

  BOOST_REQUIRE_EQUAL(cells.cells[20].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[21].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[22].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[23].size(), 1);
  BOOST_CHECK_EQUAL(cells.cells[23][0], 9);

  BOOST_REQUIRE_EQUAL(cells.cells[24].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[25].size(), 0);

  BOOST_REQUIRE_EQUAL(cells.cells[26].size(), 1);
  BOOST_CHECK_EQUAL(cells.cells[26][0], 4);
}

BOOST_AUTO_TEST_CASE(collisions_in_one_dimension) {
  const unsigned int ncell = 4;
  const double l = 20.0;
  const double lcell = l / ncell;
  const double rh2 = 1.25;
  const double rh = std::sqrt(rh2);
  NonlocalBonds nonlocal_bonds;
  UpdateConfig update_config;
  EventQueue event_queue;
  const unsigned int max_nbonds = 1;

  const std::vector<Vec3> pos = {
      // particles moving toward each other along the x-dimension
      {6.0, 10.0, 0.0},
      // particles moving away from each other along the x-dimension but will
      // collide due to periodic boundary conditions
      {2.5, 15.0, 3.0},
      // particles moving away from each other along the y-dimension; since
      // these particles are separated greater than the minimum distance l/2
      // under periodic boundary conditions, these particles do not collide
      {13.0, 11.0, 10.0},
      // particles moving toward each other along the y-dimension
      {7.0, 2.5, 6.5},
      // particles moving toward each other along the y-dimension
      {15.0, 6.0, 9.0},
      // particles moving in the same direction along the y-dimension
      {7.0, 8.0, 1.3},
      // particles moving toward each other along the y-dimension
      {1.0, 6.0, 4.2},

      {8.0, 10.0, 0.0},
      {17.5, 15.0, 3.0},
      {13.0, 18.0, 10.0},
      {7.0, 6.0, 6.5},
      {15.0, 8.0, 9.0},
      {9.5, 8.0, 1.3},
      {1.0, 8.0, 4.2}};

  const std::vector<Vec3> vel = {
      {5.0, 0.0, 0.0},  {-5.0, 0.0, 0.0}, {0.0, -8.0, 0.0}, {0.0, 2.0, 0.0},
      {0.0, 3.0, 0.0},  {-0.5, 0.0, 0.0}, {0.0, 0.0, 0.0},

      {-5.0, 0.0, 0.0}, {5.0, 0.0, 0.0},  {0.0, 3.0, 0.0},  {0.0, -6.0, 0.0},
      {0.0, -5.0, 0.0}, {-2.0, 0.0, 0.0}, {0.0, -4.0, 0.0}};

  const unsigned int nbeads = pos.size();
  std::vector<double> times(nbeads, 0);
  std::vector<uint64_t> counter(nbeads);
  std::iota(counter.begin(), counter.end(), 5);

  Cells cells(ncell, lcell);
  init_cells(pos, l, cells);

  // predict time to collision
  std::vector<double> t;
  for (unsigned int i = 0; i < 7; i++) {
    const double dvx = vel[i].x - vel[i + 7].x;
    const double dvy = vel[i].y - vel[i + 7].y;
    const double dvz = vel[i].z - vel[i + 7].z;
    const double avg_vel = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

    double dx = 0.0, dy = 0.0, dz = 0.0;
    // if particles are moving toward each other,
    if (dvx >= 0.0 && dvy >= 0.0 && dvz >= 0.0) {
      dx = pos[i].x - pos[i + 7].x;
      dy = pos[i].y - pos[i + 7].y;
      dz = pos[i].z - pos[i + 7].z;
      // else if particles are moving away from each other,
    } else if (dvx < 0.0 && dvy == 0.0 && dvz == 0.0) {
      dx = pos[i + 7].x - l - pos[i].x;
      dy = pos[i].y - pos[i + 7].y;
      dz = pos[i].z - pos[i + 7].z;
    } else if (dvx == 0.0 && dvy < 0.0 && dvz == 0.0) {
      dx = pos[i].x - pos[i + 7].x;
      dy = pos[i + 7].y - l - pos[i].y;
      dz = pos[i].z - pos[i + 7].z;
    } else if (dvx == 0.0 && dvy == 0.0 && dvz < 0.0) {
      dx = pos[i].x - pos[i + 7].x;
      dy = pos[i].y - pos[i + 7].y;
      dz = pos[i + 7].z - l - pos[i].z;
    }

    if (dx > 0.0 || dy > 0.0 || dz > 0.0) {
      std::cerr << "dx, dy, dz must be negative; check order of particles in "
                   "pos vector"
                << std::endl;
    }

    // if particles are moving toward each other along the x-dimension, ie.
    // dy and dz are 0,
    if (dx < 0.0) {
      dx += rh;
      // else if particles are moving toward each other along the
      // y-dimension, ie. dx and dz are 0,
    } else if (dy < 0.0) {
      dy += rh;
      // else if particles are moving toward each other along the
      // z-dimension, ie. dx and dy are 0,
    } else if (dz < 0.0) {
      dz += rh;
    }
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (dist < l / 2) {
      std::cout << "time " << dist / avg_vel << " of collision between " << i
                << " and " << i + 7 << std::endl;
      t.emplace_back(dist / avg_vel);
    }
  }

  add_events_for_all_beads(pos, vel, nbeads, rh2, {}, {}, l, counter,
                           event_queue, times, cells, nonlocal_bonds, {},
                           update_config, max_nbonds);

  BOOST_REQUIRE_EQUAL(event_queue.size(), 12);
  // test collision events for each particle, ie. one event occurs for a
  // collision between 0 and 7, and one for 7 and 0, but because of swap in
  // if_coll, the second event is queued as 0 colliding with 7
  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[0], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 0);
    BOOST_CHECK_EQUAL(ev.j, 7);
    BOOST_CHECK_EQUAL(ev.ni, 5);
    BOOST_CHECK_EQUAL(ev.nj, 12);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[0], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 0);
    BOOST_CHECK_EQUAL(ev.j, 7);
    BOOST_CHECK_EQUAL(ev.ni, 5);
    BOOST_CHECK_EQUAL(ev.nj, 12);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[3], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 4);
    BOOST_CHECK_EQUAL(ev.j, 11);
    BOOST_CHECK_EQUAL(ev.ni, 9);
    BOOST_CHECK_EQUAL(ev.nj, 16);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[3], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 4);
    BOOST_CHECK_EQUAL(ev.j, 11);
    BOOST_CHECK_EQUAL(ev.ni, 9);
    BOOST_CHECK_EQUAL(ev.nj, 16);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[5], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 6);
    BOOST_CHECK_EQUAL(ev.j, 13);
    BOOST_CHECK_EQUAL(ev.ni, 11);
    BOOST_CHECK_EQUAL(ev.nj, 18);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[5], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 6);
    BOOST_CHECK_EQUAL(ev.j, 13);
    BOOST_CHECK_EQUAL(ev.ni, 11);
    BOOST_CHECK_EQUAL(ev.nj, 18);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[2], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 3);
    BOOST_CHECK_EQUAL(ev.j, 10);
    BOOST_CHECK_EQUAL(ev.ni, 8);
    BOOST_CHECK_EQUAL(ev.nj, 15);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[2], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 3);
    BOOST_CHECK_EQUAL(ev.j, 10);
    BOOST_CHECK_EQUAL(ev.ni, 8);
    BOOST_CHECK_EQUAL(ev.nj, 15);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[1], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 1);
    BOOST_CHECK_EQUAL(ev.j, 8);
    BOOST_CHECK_EQUAL(ev.ni, 6);
    BOOST_CHECK_EQUAL(ev.nj, 13);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[1], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 1);
    BOOST_CHECK_EQUAL(ev.j, 8);
    BOOST_CHECK_EQUAL(ev.ni, 6);
    BOOST_CHECK_EQUAL(ev.nj, 13);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[4], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 5);
    BOOST_CHECK_EQUAL(ev.j, 12);
    BOOST_CHECK_EQUAL(ev.ni, 10);
    BOOST_CHECK_EQUAL(ev.nj, 17);
  }
  event_queue.pop();

  {
    const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, t[4], 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 5);
    BOOST_CHECK_EQUAL(ev.j, 12);
    BOOST_CHECK_EQUAL(ev.ni, 10);
    BOOST_CHECK_EQUAL(ev.nj, 17);
  }
  event_queue.pop();
}

BOOST_AUTO_TEST_CASE(cell_crossing_time) {
  const unsigned int ncell = 3;
  const double lcell = 1.0;
  Cells cells(ncell, lcell);
  BeadCellEvent::Wall wall;
  unsigned int ixn, iyn, izn;

  // vxi = 0, vyi = 0, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(!t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 0.0, 0.0, ixo,
                              iyo, izo, dt, wall, ixn, iyn, izn));
  }
  // vxi > 0, vyi = 0, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 5.0, 0.0, 0.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi = 0, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -5.0, 0.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xneg);
    BOOST_CHECK_EQUAL(ixn, 1);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi = 0, vyi > 0, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 5.0, 0.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.2);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::ypos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 0);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi = 0, vyi < 0 (at the boundary), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.0, 3.0, cells, 0.0, -5.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi = 0, vyi < 0 (in the center), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.5, 3.0, cells, 0.0, -5.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi = 0, vyi = 0, vzi > 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 0.0, 5.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zpos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 0);
  }
  // vxi = 0, vyi = 0, vzi < 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 0.0, -5.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 1);
  }
  // vxi > 0, vyi > 0, vxi > vyi, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 5.0, 4.0, 0.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi > 0, vyi > 0, vxi < vyi, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 4.0, 5.0, 0.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.125);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi > 0, vyi < 0, vxi > vyi (at the boundary), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 5.0, -4.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi > 0, vyi < 0, vxi > vyi (in the center), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.5, 3.0, cells, 5.0, -4.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi > 0, vxi < vyi, vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -4.0, 5.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.125);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xneg);
    BOOST_CHECK_EQUAL(ixn, 1);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi > 0, vzi > 0, vxi > vzi, vyi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 5.0, 0.0, 4.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi > 0, vzi > 0, vxi < vzi, vyi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 4.0, 0.0, 5.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zpos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 0);
  }
  // vxi > 0, vzi < 0, vxi > vzi, vyi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 5.0, 0.0, -4.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xpos);
    BOOST_CHECK_EQUAL(ixn, 0);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vzi > 0, vxi < vzi, vyi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -4.0, 0.0, 5.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zpos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 0);
  }
  // vyi > 0, vzi > 0, vyi > vzi, vxi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 5.0, 4.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.125);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zpos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 0);
  }
  // vyi > 0, vzi > 0, vyi < vzi, vxi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 4.0, 5.0, ixo, iyo,
                             izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zpos);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 0);
  }
  // vyi > 0, vzi < 0, vyi > vzi, vxi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, 5.0, -4.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.125);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 1);
  }
  // vyi < 0, vzi > 0, vyi < vzi, vxi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, 0.0, -4.0, 5.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vxi > vyi (at the boundary), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -4.0, -5.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vxi > vyi (in the center), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.0, 3.0, cells, -4.0, -5.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vzi < 0, vxi > vyi (at the boundary) > vzi
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -4.0, -5.0, -6.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vzi < 0, vxi > vyi (in the center) > vzi
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.0, 3.0, cells, -4.0, -5.0, -6.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::zneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 1);
  }
  // vxi < 0, vyi < 0, vxi < vyi (at the boundary), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -5.0, -4.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vxi < vyi (in the center), vzi = 0
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.5, 3.0, cells, -5.0, -4.0, 0.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xneg);
    BOOST_CHECK_EQUAL(ixn, 1);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vzi < 0, vxi < vyi (at the boundary) < vzi
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.0, 2.5, 3.0, cells, -5.0, -4.0, -3.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.0);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::yneg);
    BOOST_CHECK_EQUAL(ixn, 2);
    BOOST_CHECK_EQUAL(iyn, 1);
    BOOST_CHECK_EQUAL(izn, 2);
  }
  // vxi < 0, vyi < 0, vzi < 0, vxi < vyi (in the center) < vzi
  {
    double dt = 0.0;
    unsigned int ixo = 2;
    unsigned int iyo = 2;
    unsigned int izo = 2;
    BOOST_CHECK(t_until_cell(2.5, 2.5, 2.5, 3.0, cells, -5.0, -4.0, -3.0, ixo,
                             iyo, izo, dt, wall, ixn, iyn, izn));
    BOOST_CHECK_EQUAL(dt, 0.1);
    BOOST_CHECK_EQUAL(wall, BeadCellEvent::xneg);
    BOOST_CHECK_EQUAL(ixn, 1);
    BOOST_CHECK_EQUAL(iyn, 2);
    BOOST_CHECK_EQUAL(izn, 2);
  }
}

BOOST_AUTO_TEST_CASE(cell_position_with_vel) {
  const unsigned int ncell = 3;
  const double lcell = 1.0;

  const std::vector<Vec3> pos = {
      {1.8, 0.5, 2.3}, {0.0, 0.5, 0.0}, {2.0, 2.3, 1.4}, {0.5, 0.5, 2.7},
      {2.0, 2.0, 2.0}, {0.5, 1.8, 0.2}, {2.2, 2.6, 1.9}, {1.5, 2.0, 3.0},
      {1.0, 2.7, 0.6}, {2.5, 1.5, 2.5}};

  const std::vector<Vec3> vel = {
      {5.0, 0.0, 0.0},  {0.0, -5.0, 0.0}, {0.0, 5.0, 0.0},  {3.0, -2.0, 0.0},
      {4.0, 0.0, 0.0},  {-4.0, 0.0, 0.0}, {1.0, -6.0, 0.0}, {2.0, 1.0, 0.0},
      {-2.0, 2.0, 0.0}, {4.0, -3.0, 0.0}};

  Cells cells(ncell, lcell);
  init_cells(pos, ncell, cells);

  BOOST_REQUIRE_EQUAL(cells.cell_x.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_y.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_z.size(), 10);
  BOOST_CHECK_EQUAL(cells.cell_x[0], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[0], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[0], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[2], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[2], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[2], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[3], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[3], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[3], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[5], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[6], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[6], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[6], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[7], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[7], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[7], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[8], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[8], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[8], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[9], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[9], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[9], 2);

  {
    BOOST_REQUIRE_EQUAL(cells.cells[19].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[19][0], 0);
    BOOST_REQUIRE_EQUAL(cells.cells[20].size(), 0);

    move_to_new_cell(cells, 0, 2, 0, 2);

    BOOST_REQUIRE_EQUAL(cells.cells[20].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[20][0], 0);
    BOOST_REQUIRE_EQUAL(cells.cells[19].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[0][0], 1);
    BOOST_REQUIRE_EQUAL(cells.cells[6].size(), 0);

    move_to_new_cell(cells, 1, 0, 2, 0);

    BOOST_REQUIRE_EQUAL(cells.cells[6].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[6][0], 1);
    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[17].size(), 2);
    BOOST_CHECK_EQUAL(cells.cells[17][0], 2);
    BOOST_CHECK_EQUAL(cells.cells[17][1], 6);
    BOOST_REQUIRE_EQUAL(cells.cells[11].size(), 0);

    move_to_new_cell(cells, 2, 2, 0, 1);

    BOOST_REQUIRE_EQUAL(cells.cells[11].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[11][0], 2);
    BOOST_REQUIRE_EQUAL(cells.cells[17].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[17][0], 6);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[18].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[18][0], 3);
    BOOST_REQUIRE_EQUAL(cells.cells[19].size(), 0);

    move_to_new_cell(cells, 3, 1, 0, 2);

    BOOST_REQUIRE_EQUAL(cells.cells[19].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[19][0], 3);
    BOOST_REQUIRE_EQUAL(cells.cells[18].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[26].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[26][0], 4);
    BOOST_REQUIRE_EQUAL(cells.cells[24].size(), 0);

    move_to_new_cell(cells, 4, 0, 2, 2);

    BOOST_REQUIRE_EQUAL(cells.cells[24].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[24][0], 4);
    BOOST_REQUIRE_EQUAL(cells.cells[26].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[3].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[3][0], 5);
    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 0);

    move_to_new_cell(cells, 5, 0, 0, 0);

    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[0][0], 5);
    BOOST_REQUIRE_EQUAL(cells.cells[3].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[17].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[17][0], 6);
    BOOST_REQUIRE_EQUAL(cells.cells[14].size(), 0);

    move_to_new_cell(cells, 6, 2, 1, 1);

    BOOST_REQUIRE_EQUAL(cells.cells[14].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[14][0], 6);
    BOOST_REQUIRE_EQUAL(cells.cells[17].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[7].size(), 2);
    BOOST_CHECK_EQUAL(cells.cells[7][0], 7);
    BOOST_CHECK_EQUAL(cells.cells[7][1], 8);
    BOOST_REQUIRE_EQUAL(cells.cells[8].size(), 0);

    move_to_new_cell(cells, 7, 2, 2, 0);

    BOOST_REQUIRE_EQUAL(cells.cells[8].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[8][0], 7);
    BOOST_REQUIRE_EQUAL(cells.cells[7].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[7][0], 8);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[7].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[7][0], 8);
    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[0][0], 5);

    move_to_new_cell(cells, 8, 0, 0, 0);

    BOOST_REQUIRE_EQUAL(cells.cells[0].size(), 2);
    BOOST_CHECK_EQUAL(cells.cells[0][0], 5);
    BOOST_CHECK_EQUAL(cells.cells[0][1], 8);
    BOOST_REQUIRE_EQUAL(cells.cells[7].size(), 0);
  }
  {
    BOOST_REQUIRE_EQUAL(cells.cells[23].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[23][0], 9);
    BOOST_REQUIRE_EQUAL(cells.cells[21].size(), 0);

    move_to_new_cell(cells, 9, 0, 1, 2);

    BOOST_REQUIRE_EQUAL(cells.cells[21].size(), 1);
    BOOST_CHECK_EQUAL(cells.cells[21][0], 9);
    BOOST_REQUIRE_EQUAL(cells.cells[23].size(), 0);
  }

  BOOST_REQUIRE_EQUAL(cells.cell_x.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_y.size(), 10);
  BOOST_REQUIRE_EQUAL(cells.cell_z.size(), 10);
  BOOST_CHECK_EQUAL(cells.cell_x[0], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[0], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[0], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[1], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[1], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[2], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[2], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[2], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[3], 1);
  BOOST_CHECK_EQUAL(cells.cell_y[3], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[3], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[4], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[4], 2);
  BOOST_CHECK_EQUAL(cells.cell_x[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[5], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[6], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[6], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[6], 1);
  BOOST_CHECK_EQUAL(cells.cell_x[7], 2);
  BOOST_CHECK_EQUAL(cells.cell_y[7], 2);
  BOOST_CHECK_EQUAL(cells.cell_z[7], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[8], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[8], 0);
  BOOST_CHECK_EQUAL(cells.cell_z[8], 0);
  BOOST_CHECK_EQUAL(cells.cell_x[9], 0);
  BOOST_CHECK_EQUAL(cells.cell_y[9], 1);
  BOOST_CHECK_EQUAL(cells.cell_z[9], 2);
}

BOOST_AUTO_TEST_CASE(test_add_events_for_bead_after_crossing) {
  const unsigned int ncell = 3;
  const double l = 6.0;
  const double lcell = l / ncell;
  const double rh2 = 1.0;
  NonlocalBonds nonlocal_bonds;
  UpdateConfig update_config;
  const unsigned int max_nbonds = 1;

  std::vector<Vec3> pos{
      {3.0, 3.0, 3.0},
      // nearest neighbor
      {3.0, 2.0, 3.0},
      // next-nearest neighbor
      {3.0, 4.0, 3.0},

      // particles approach i = 0 from left side
      // z = 0
      {1.0, 1.0, 1.0},
      {1.0, 3.0, 1.0},
      {1.0, 5.0, 1.0},
      // z = 1
      {1.0, 1.0, 3.0},
      {1.0, 3.0, 3.0},
      {1.0, 5.0, 3.0},
      // z = 2
      {1.0, 1.0, 5.0},
      {1.0, 3.0, 5.0},
      {1.0, 5.0, 5.0},

      // particles approach i = 0 from right side
      // z = 0
      {5.0, 1.0, 1.0},
      {5.0, 3.0, 1.0},
      {5.0, 5.0, 1.0},
      // z = 1
      {5.0, 1.0, 3.0},
      {5.0, 3.0, 3.0},
      {5.0, 5.0, 3.0},
      // z = 2
      {5.0, 1.0, 5.0},
      {5.0, 3.0, 5.0},
      {5.0, 5.0, 5.0},

      // particles approach i = 0 from below
      // z = 0
      {3.0, 1.0, 1.0},
      // z = 1
      {3.0, 1.0, 3.0},
      // z = 2
      {3.0, 1.0, 5.0},

      // particles approach i = 0 from above
      // z = 0
      {3.0, 5.0, 1.0},
      // z = 1
      {3.0, 5.0, 3.0},
      // z = 2
      {3.0, 5.0, 5.0},

      // particles approach i = 0 from the back
      // z = 0
      {3.0, 3.0, 1.0},

      // particles approach i = 0 from the front
      // z = 2
      {3.0, 3.0, 5.0},
  };

  std::vector<Vec3> vel{
      {0.0, 0.0, 0.0},
      // nearest neighbor
      {-3.0, -3.0, 0.0},
      // next-nearest neighbor
      {-4.0, 4.0, 0.0},
      // particles approach i = 0 from left side
      // z = 0
      {5.0, 5.0, 5.0},
      {6.0, 0.0, 6.0},
      {2.0, -2.0, 2.0},
      // z = 1
      {3.0, 3.0, 0.0},
      {1.0, 0.0, 0.0},
      {4.0, -4.0, 0.0},
      // z = 2
      {8.0, 8.0, -8.0},
      {2.0, 0.0, -2.0},
      {7.0, -7.0, -7.0},

      // particles approach i = 0 from right side
      // z = 0
      {-4.0, 4.0, 4.0},
      {-5.0, 0.0, 5.0},
      {-6.0, -6.0, 6.0},
      // z = 1
      {-7.0, 7.0, 0.0},
      {-2.0, 0.0, 0.0},
      {-1.0, -1.0, 0.0},
      // z = 2
      {-9.0, 9.0, -9.0},
      {-8.0, 0.0, -8.0},
      {-3.0, -3.0, -3.0},

      // particles approach i = 0 from below
      // z = 0
      {0.0, 5.5, 5.5},
      // z = 1
      {0.0, 8.0, 0.0},
      // z = 2
      {0.0, 1.5, -1.5},

      // particles approach i = 0 from above
      // z = 0
      {0.0, -9.0, 9.0},
      // z = 1
      {0.0, -2.5, 0.0},
      // z = 2
      {0.0, -3.5, -3.5},

      // particles approach i = 0 from the back
      // z = 0
      {0.0, 0.0, 10.0},

      // particles approach i = 0 from the front
      // z = 2
      {0.0, 0.0, -5.0},
  };

  const unsigned int nbeads = pos.size();
  std::vector<double> times(nbeads, 0);
  std::vector<uint64_t> counter(nbeads);
  std::iota(counter.begin(), counter.end(), 5);

  Cells cells(ncell, lcell);
  init_cells(pos, l, cells);

  // particles approach i = 0 from left side
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::xneg, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.17783121635, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 9);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 14);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.20323567583, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 11);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 16);
    }
    event_queue.pop();

    // to find expected d value, treat each particle as a point particle
    // and find the distance between them, and then subtract the radii of
    // each particle from this distance
    // t = d / v = (sqrt(8) - 0.5 - 0.5) / (6 * sqrt(2))
    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.21548220313, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 4);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 9);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.28452994616, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 3);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 8);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.3232233047, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 8);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 13);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.43096440627, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 6);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 11);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.6464466094, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 10);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 15);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.7113248654, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 5);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 10);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 1.0, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 7);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 12);
    }
    event_queue.pop();
  }

  // particles approach i = 0 from right side
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::xpos, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.15807219231, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 18);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 23);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.16161165235, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 19);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 24);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.18469903125, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 15);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 20);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.23710828846, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 14);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 19);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.25857864376, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 13);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 18);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.3556624327, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 12);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 17);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.47421657693, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 20);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 25);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.5, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 16);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 21);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 1.29289321881, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 17);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 22);
    }
    event_queue.pop();
  }

  // particles approach i = 0 from below
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::yneg, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.125, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 22);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 27);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.15807219231, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 18);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 23);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.17783121635, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 9);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 14);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.18469903125, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 15);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 20);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.23507149433, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 21);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 26);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.28452994616, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 3);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 8);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.3556624327, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 12);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 17);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.43096440627, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 6);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 11);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.86192881254, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 23);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 28);
    }
    event_queue.pop();
  }

  // particles approach i = 0 from above
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::ypos, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.14365480209, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 24);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 29);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.20323567583, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 11);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 16);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.23710828846, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 14);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 19);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.3232233047, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 8);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 13);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.36939806251, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 26);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 31);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.4, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 25);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 30);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.47421657693, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 20);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 25);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.7113248654, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 5);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 10);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 1.29289321881, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 17);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 22);
    }
    event_queue.pop();
  }

  // particles approach i = 0 from the back
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::zneg, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.1, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 27);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 32);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.14365480209, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 24);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 29);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.21548220313, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 4);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 9);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.23507149433, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 21);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 26);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.23710828846, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 14);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 19);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.25857864376, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 13);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 18);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.28452994616, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 3);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 8);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.3556624327, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 12);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 17);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.7113248654, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 5);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 10);
    }
    event_queue.pop();
  }

  // particles approach i = 0 from the front
  {
    EventQueue event_queue;
    add_events_for_bead_after_crossing(
        pos, vel, rh2, {}, {}, l, counter, event_queue, times, cells, 0,
        BeadCellEvent::zpos, nonlocal_bonds, {}, update_config, max_nbonds);
    BOOST_REQUIRE_EQUAL(event_queue.size(), 9);

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.15807219231, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 18);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 23);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.16161165235, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 19);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 24);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.17783121635, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 9);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 14);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.2, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 28);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 33);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.20323567583, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 11);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 16);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.36939806251, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 26);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 31);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.47421657693, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 20);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 25);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.6464466094, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 10);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 15);
    }
    event_queue.pop();

    {
      const MinNonlocalInnerEvent &ev =
          std::get<MinNonlocalInnerEvent>(event_queue.top());
      BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.86192881254, 1e-10);
      BOOST_CHECK_EQUAL(ev.i, 0);
      BOOST_CHECK_EQUAL(ev.j, 23);
      BOOST_CHECK_EQUAL(ev.ni, 5);
      BOOST_CHECK_EQUAL(ev.nj, 28);
    }
    event_queue.pop();
  }
}

BOOST_AUTO_TEST_CASE(points_on_unit_sphere) {
  int seed = 40;
  std::random_device rd;
  std::mt19937 mt(seed);

  {
    for (unsigned int i = 0; i < 1000000; i++) {
      double x;
      double y;
      double z;
      unit_sphere(mt, x, y, z);
      double dist = std::sqrt(x * x + y * y + z * z);
      BOOST_CHECK_CLOSE_FRACTION(dist, 1, 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(initialize_bonds) {
  const double l = 3.0;
  const double near_min2 = 1.0 * 1.0;
  const double near_max2 = 1.17 * 1.17;

  const std::vector<Vec3> pos = {{2.0, 1.0, 0.0},    {1.0, 1.0, 0.0},
                                 {1.0, 1.0, 1.17},   {2.17, 1.0, 1.17},
                                 {2.17, 2.0, 1.17},  {2.17, 0.17, 1.17},
                                 {0.17, 0.17, 1.17}, {0.17, 0.17, 2.17}};

  const std::vector<Vec3> vel = {
      {5.0, 0.0, 0.0}, {-5.0, 0.0, 0.0}, {0.0, 0.0, 3.0},  {-5.0, 0.0, 0.0},
      {0.0, 5.0, 0.0}, {-5.0, 0.0, 0.0}, {0.0, 0.0, -1.0}, {2.0, 0.0, 0.0}};

  const unsigned int nbeads = pos.size();
  std::vector<double> times(nbeads, 0);
  std::vector<uint64_t> counter(nbeads);
  std::iota(counter.begin(), counter.end(), 5);

  EventQueue event_queue;
  init_nearest_bond_events(pos, vel, nbeads, l, counter, event_queue, times,
                           near_min2, near_max2);

  BOOST_REQUIRE_EQUAL(event_queue.size(), 7);
  {
    const MaxNearestEvent &ev = std::get<MaxNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.0, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 1);
    BOOST_CHECK_EQUAL(ev.j, 2);
    BOOST_CHECK_EQUAL(ev.ni, 6);
    BOOST_CHECK_EQUAL(ev.nj, 7);
  }
  event_queue.pop();

  {
    const MaxNearestEvent &ev = std::get<MaxNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.017, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 0);
    BOOST_CHECK_EQUAL(ev.j, 1);
    BOOST_CHECK_EQUAL(ev.ni, 5);
    BOOST_CHECK_EQUAL(ev.nj, 6);
  }
  event_queue.pop();

  {
    const MaxNearestEvent &ev = std::get<MaxNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.031825642422102, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 3);
    BOOST_CHECK_EQUAL(ev.j, 4);
    BOOST_CHECK_EQUAL(ev.ni, 8);
    BOOST_CHECK_EQUAL(ev.nj, 9);
  }
  event_queue.pop();

  {
    const MaxNearestEvent &ev = std::get<MaxNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.033901746176141, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 5);
    BOOST_CHECK_EQUAL(ev.j, 6);
    BOOST_CHECK_EQUAL(ev.ni, 10);
    BOOST_CHECK_EQUAL(ev.nj, 11);
  }
  event_queue.pop();

  {
    const MinNearestEvent &ev = std::get<MinNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.035112707285374, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 2);
    BOOST_CHECK_EQUAL(ev.j, 3);
    BOOST_CHECK_EQUAL(ev.ni, 7);
    BOOST_CHECK_EQUAL(ev.nj, 8);
  }
  event_queue.pop();

  {
    const MinNearestEvent &ev = std::get<MinNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.037558197402123, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 4);
    BOOST_CHECK_EQUAL(ev.j, 5);
    BOOST_CHECK_EQUAL(ev.ni, 9);
    BOOST_CHECK_EQUAL(ev.nj, 10);
  }
  event_queue.pop();

  {
    const MaxNearestEvent &ev = std::get<MaxNearestEvent>(event_queue.top());
    BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.137312911107772, 1e-10);
    BOOST_CHECK_EQUAL(ev.i, 6);
    BOOST_CHECK_EQUAL(ev.j, 7);
    BOOST_CHECK_EQUAL(ev.ni, 11);
    BOOST_CHECK_EQUAL(ev.nj, 12);
  }
  event_queue.pop();
}

BOOST_AUTO_TEST_CASE(struct_nonlocal_bonds) {
  {
    const unsigned int n = 59;
    const NonlocalBonds nb{
        {{6, 18, 1.5}, {18, 30, 1.5}, {22, 34, 1.5}, {34, 46, 1.5}, {42, 54, 1.5}, {46, 58, 1.5}}};

    const std::vector<std::pair<unsigned int, unsigned int>> bonds = {
        {6, 18}, {18, 30}, {22, 34}, {34, 46}, {42, 54}, {46, 58}};

    // loop through all possible bond pairs possible with n indices (beads)
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = i + 1; j < n; j++) {
          // find output given by get bond mask (first element is config and second is the rc)
          const std::tuple<Config, double> bond_mask_tuple = nb.get_bond_mask(i, j);
          // if the (i, j) pair is in bonds then bond mask is true else false
          if (std::find(bonds.begin(), bonds.end(), std::make_pair(i, j)) != bonds.end()) {
              BOOST_CHECK_MESSAGE(std::get<0>(bond_mask_tuple), "check if nb(" << i << ", " << j << ") has failed");
          }
          else {
              BOOST_CHECK_MESSAGE(!std::get<0>(bond_mask_tuple), "check if nb(" << i << ", " << j << ") has failed");
          }
      }
    }
  }

  {
    const unsigned int n = 47;

    const NonlocalBonds nb{
        {{2, 10, 1.5}, {10, 18, 1.5}, {14, 22, 1.5}, {26, 34, 1.5}, {30, 38, 1.5}, {38, 46, 1.5}}
    };

    const std::vector<std::pair<unsigned int, unsigned int>> bonds = {
            {2, 10}, {10, 18}, {14, 22}, {26, 34}, {30, 38}, {38, 46}
    };

    // loop through all possible bond pairs possible with n indices (beads)
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            // find output given by get bond mask (first element is config and second is the rc)
            const std::tuple<Config, double> bond_mask_tuple = nb.get_bond_mask(i, j);
            // if the (i, j) pair is in bonds then bond mask is true else false
            if (std::find(bonds.begin(), bonds.end(), std::make_pair(i, j)) != bonds.end()) {
                BOOST_CHECK_MESSAGE(std::get<0>(bond_mask_tuple), "check nb(" << i << ", " << j << ") has failed");
            }
            else {
                BOOST_CHECK_MESSAGE(!std::get<0>(bond_mask_tuple), "check nb(" << i << ", " << j << ") has failed");
            }
        }
    }
  }
}

BOOST_AUTO_TEST_CASE(nonlocal_bonding_events) {
    const unsigned int ncell = 4;
    const double l = 20.0;
    const double lcell = l / ncell;
    const double rc = std::sqrt(3.0);
    const NonlocalBonds permanent_bonds{{{2, 6, rc}}};
    const NonlocalBonds transient_bonds{{{6, 10, rc}}};
    UpdateConfig update_config;
    const double rh2 = 1.0;
    EventQueue event_queue;
    const unsigned int max_nbonds = 1;

    const std::vector<Vec3> pos = {
            {6.0, 10.0, 0.0}, {2.5, 15.0, 3.0},  {9.0, 10.0, 0.0}, {13.0, 11.0, 10.0},
            {7.0, 6.0, 6.5},  {15.0, 8.0, 9.0},  {8.0, 10.0, 0.0}, {9.5, 8.0, 1.3},
            {1.0, 8.0, 4.2},  {17.5, 15.0, 3.0}, {9.5, 10.0, 0.0}};

    const std::vector<Vec3> vel = {
            {5.0, 0.0, 0.0},  {-5.0, 0.0, 0.0}, {5.0, 0.0, 0.0},  {0.0, -8.0, 0.0},
            {0.0, -6.0, 0.0}, {0.0, -5.0, 0.0}, {-5.0, 0.0, 0.0}, {-2.0, 0.0, 0.0},
            {0.0, -4.0, 0.0}, {5.0, 0.0, 0.0},  {5.0, 0.0, 0.0}};

    const unsigned int nbeads = pos.size();
    std::vector<double> times(nbeads, 0);
    std::vector<uint64_t> counter(nbeads);
    std::iota(counter.begin(), counter.end(), 5);

    Cells cells(ncell, lcell);
    init_cells(pos, l, cells);

    const Config p_bond_mask = std::get<0>(permanent_bonds.get_bond_mask(2, 6));
    update_config.flip_bond(p_bond_mask);
    BOOST_CHECK(update_config.bonded(p_bond_mask));
    const Config t_bond_mask = std::get<0>(transient_bonds.get_bond_mask(6, 10));
    BOOST_CHECK(update_config.bonded(t_bond_mask));

for (unsigned int i = 0; i < nbeads; i++) {
    for (unsigned int j = 0; j < nbeads; j++) {
        if_coll(pos, vel, rh2, {}, {}, l, counter, event_queue, times, i, j,
                transient_bonds, permanent_bonds, update_config, max_nbonds);
    }
}

BOOST_REQUIRE_EQUAL(event_queue.size(), 8);
{
const MaxNonlocalOuterEvent &ev = std::get<MaxNonlocalOuterEvent>(event_queue.top());
std::cout << ev << std::endl;
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.0232050807568877, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 6);
BOOST_CHECK_EQUAL(ev.j, 10);
BOOST_CHECK_EQUAL(ev.ni, 11);
BOOST_CHECK_EQUAL(ev.nj, 15);
}
event_queue.pop();

{
const MaxNonlocalOuterEvent &ev =
        std::get<MaxNonlocalOuterEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.0232050807568877, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 6);
BOOST_CHECK_EQUAL(ev.j, 10);
BOOST_CHECK_EQUAL(ev.ni, 11);
BOOST_CHECK_EQUAL(ev.nj, 15);
}
event_queue.pop();

{
const MaxNonlocalOuterEvent &ev =
        std::get<MaxNonlocalOuterEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.0732050807568877, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 2);
BOOST_CHECK_EQUAL(ev.j, 6);
BOOST_CHECK_EQUAL(ev.ni, 7);
BOOST_CHECK_EQUAL(ev.nj, 11);
}
event_queue.pop();

{
const MaxNonlocalOuterEvent &ev =
        std::get<MaxNonlocalOuterEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.0732050807568877, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 2);
BOOST_CHECK_EQUAL(ev.j, 6);
BOOST_CHECK_EQUAL(ev.ni, 7);
BOOST_CHECK_EQUAL(ev.nj, 11);
}
event_queue.pop();

{
const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.1, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 0);
BOOST_CHECK_EQUAL(ev.j, 6);
BOOST_CHECK_EQUAL(ev.ni, 5);
BOOST_CHECK_EQUAL(ev.nj, 11);
}
event_queue.pop();

{
const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.1, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 0);
BOOST_CHECK_EQUAL(ev.j, 6);
BOOST_CHECK_EQUAL(ev.ni, 5);
BOOST_CHECK_EQUAL(ev.nj, 11);
}
event_queue.pop();

{
const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.4, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 1);
BOOST_CHECK_EQUAL(ev.j, 9);
BOOST_CHECK_EQUAL(ev.ni, 6);
BOOST_CHECK_EQUAL(ev.nj, 14);
}
event_queue.pop();

{
const MinNonlocalInnerEvent &ev =
        std::get<MinNonlocalInnerEvent>(event_queue.top());
BOOST_CHECK_CLOSE_FRACTION(ev.t, 0.4, 1e-10);
BOOST_CHECK_EQUAL(ev.i, 1);
BOOST_CHECK_EQUAL(ev.j, 9);
BOOST_CHECK_EQUAL(ev.ni, 6);
BOOST_CHECK_EQUAL(ev.nj, 14);
}
event_queue.pop();
}

BOOST_AUTO_TEST_CASE(test_rodrigues) {
  {
    const std::vector<Vec3> pos = {
        {1.0, 2.0, 0.0}, {2.0, 2.0, 0.0}, {3.0, 2.5, 0.0}, {4.0, 2.5, 0.0}};

    const unsigned int nbeads = pos.size();
    const double theta = M_PI;
    const unsigned int ind = 1;
    const double l = 5.0;

    std::vector<Vec3> pos_trial;
    for (unsigned int i = 0; i < nbeads; i++) {
      pos_trial.emplace_back(pos[i]);
    }

    rodrigues_rotation(pos, theta, pos_trial, ind, l);

    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].x, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].z, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].x, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].z, 0.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[2].x, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[2].y, 1.5, 1e-10);
    BOOST_CHECK_SMALL(pos_trial[2].z, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[3].x, 4.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[3].y, 1.5, 1e-10);
    BOOST_CHECK_SMALL(pos_trial[3].z, 1e-10);
  }

  {
    const std::vector<Vec3> pos = {{1.0, 2.0, 1.0}, {2.0, 2.0, 1.0},
                                   {2.0, 2.0, 0.0}, {2.0, 2.0, 3.0},
                                   {3.0, 2.5, 1.0}, {3.0, 1.5, 1.0}};

    const unsigned int nbeads = pos.size();
    const double theta = M_PI / 2;
    const unsigned int ind = 1;
    const double l = 4.0;

    std::vector<Vec3> pos_trial;
    for (unsigned int i = 0; i < nbeads; i++) {
      pos_trial.emplace_back(pos[i]);
    }

    rodrigues_rotation(pos, theta, pos_trial, ind, l);

    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].x, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[0].z, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].x, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[1].z, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[2].x, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[2].y, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[2].z, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[3].x, 2.0, 1e-10);
    BOOST_CHECK_SMALL(pos_trial[3].y, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[3].z, 1.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[4].x, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[4].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[4].z, 1.5, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[5].x, 3.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[5].y, 2.0, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(pos_trial[5].z, 0.5, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(test_trial_config) {
  const double l = 4.0;
  const NonlocalBonds nonlocal_bonds{{{2, 6, 1.5}, {6, 10, 1.5}}};

  int seed = 40;
  std::random_device rd;
  std::mt19937 mt(seed);

  std::vector<Vec3> pos = {{0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0},
                           {1.0, 2.0, 1.0}, {1.0, 2.0, 2.0}, {2.0, 1.0, 2.0},
                           {2.0, 1.0, 1.0}, {2.0, 2.0, 1.0}, {2.0, 2.0, 0.0},
                           {2.0, 1.0, 0.0}, {3.0, 1.0, 1.0}};

  const unsigned int nbeads = pos.size();
  double theta = M_PI / 2;
  unsigned int ind = 6;

  std::vector<double> s_bias;
  const unsigned int nbonds = nonlocal_bonds.get_nbonds();
  init_s(s_bias, nbonds);

  UpdateConfig orig_config = config_int(pos, l, nonlocal_bonds);
  BOOST_CHECK_EQUAL(orig_config.config, 3);

  std::vector<Vec3> pos_trial;
  for (unsigned int i = 0; i < nbeads; i++) {
    pos_trial.emplace_back(pos[i]);
  }

  rodrigues_rotation(pos, theta, pos_trial, ind, l);

  UpdateConfig trial_config = config_int(pos_trial, l, nonlocal_bonds);

  if (accept_move(s_bias, orig_config, trial_config, mt)) {
    std::swap(pos, pos_trial);
    std::swap(orig_config, trial_config);
  }

  BOOST_CHECK_EQUAL(trial_config.config, 3);
  BOOST_CHECK(accept_move(s_bias, orig_config, trial_config, mt));
}

BOOST_AUTO_TEST_CASE(test_compute_entropy) {
  {
    // each state is visited at least once
    std::vector<uint64_t> config_int;
    config_int = {3, 5, 0, 1, 6, 5, 2, 1, 3, 4, 2, 1, 4, 1, 4, 7, 6, 0, 5, 4};

    const unsigned int nbonds = 3;
    const unsigned int nstates = pow(2, nbonds);
    std::vector<uint64_t> config_count(nstates);
    unsigned int config_int_size = config_int.size();
    const double pos_scale = 1.0;
    const double neg_scale = 2.0;

    for (unsigned int i = 0; i < config_int_size; i++) {
      config_count[config_int[i]]++;
    }

    std::vector<double> s_bias;
    init_s(s_bias, nbonds);

    compute_entropy(config_count, s_bias, nstates, pos_scale, neg_scale);

    BOOST_CHECK_CLOSE_FRACTION(s_bias[0], 9.693147180559945, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[1], 7.38629436111989, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[2], 6.693147180559945, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[3], 3.693147180559945, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[4], 7.38629436111989, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[5], 4.09861228866811, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[6], 3.693147180559945, 1e-10);
    BOOST_CHECK_SMALL(s_bias[7], 1e-10);
  }

  {
    // each state is visited at least once except the native state
    std::vector<uint64_t> config_int;
    config_int = {3, 5, 0, 0, 6, 5, 2, 1, 3, 4};

    const unsigned int nbonds = 3;
    const unsigned int nstates = pow(2, nbonds);
    std::vector<uint64_t> config_count(nstates);
    unsigned int config_int_size = config_int.size();
    const double pos_scale = 1.0;
    const double neg_scale = 2.0;

    for (unsigned int i = 0; i < config_int_size; i++) {
      config_count[config_int[i]]++;
    }

    std::vector<double> s_bias;
    init_s(s_bias, nbonds);

    compute_entropy(config_count, s_bias, nstates, pos_scale, neg_scale);

    BOOST_CHECK_CLOSE_FRACTION(s_bias[0], 9.916290731874155, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[1], 6.22314355131421, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[2], 6.22314355131421, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[3], 3.916290731874155, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[4], 6.22314355131421, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[5], 3.916290731874155, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[6], 3.22314355131421, 1e-10);
    BOOST_CHECK_SMALL(s_bias[7], 1e-10);

    config_int.clear();
    config_count.clear();
    config_count.resize(nstates);

    // each state is visited at least once except the first state
    config_int = {7, 4, 4, 6, 6, 3, 5, 2, 3, 1};

    for (unsigned int i = 0; i < config_int_size; i++) {
      config_count[config_int[i]]++;
    }

    compute_entropy(config_count, s_bias, nstates, pos_scale, neg_scale);

    BOOST_CHECK_CLOSE_FRACTION(s_bias[0], 9.693147180559946, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[1], 6.22314355131421, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[2], 6.22314355131421, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[3], 4.6094379124341, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[4], 6.916290731874155, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[5], 3.916290731874155, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[6], 3.916290731874155, 1e-10);
    BOOST_CHECK_SMALL(s_bias[7], 1e-10);
  }

  {
    // some states are visited, but several are not
    std::vector<uint64_t> config_int;
    config_int = {7, 3, 3, 7, 0, 6, 6, 6, 1, 3};

    const unsigned int nbonds = 3;
    const unsigned int nstates = pow(2, nbonds);
    std::vector<uint64_t> config_count(nstates);
    unsigned int config_int_size = config_int.size();
    const double pos_scale = 1.0;
    const double neg_scale = 2.0;

    for (unsigned int i = 0; i < config_int_size; i++) {
      config_count[config_int[i]]++;
    }

    std::vector<double> s_bias;
    init_s(s_bias, nbonds);

    compute_entropy(config_count, s_bias, nstates, pos_scale, neg_scale);

    BOOST_CHECK_CLOSE_FRACTION(s_bias[0], 8.30685281944005, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[1], 5.30685281944005, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[2], 5.08370926812585, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[3], 3.405465108108164, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[4], 5.08370926812585, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[5], 2.08370926812584, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(s_bias[6], 3.405465108108164, 1e-10);
    BOOST_CHECK_SMALL(s_bias[7], 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(test_g_val) {
  std::vector<uint64_t> config_count;
  config_count = {105, 95, 102, 98, 98, 102};

  const unsigned int nstates = config_count.size();

  double g_val = compute_g_val(config_count, nstates);

  BOOST_CHECK_CLOSE_FRACTION(g_val, 0.660219210319822, 1e-10);
}

BOOST_AUTO_TEST_CASE(check_dist_constraints) {
  const double near_min2 = 1.0 * 1.0;
  const double near_max2 = 1.17 * 1.17;
  const double nnear_min2 = 1.4 * 1.4;
  const double nnear_max2 = 1.67 * 1.67;
  const double l = 5.0;

  const std::vector<Vec3> pos_0 = {
      {1.5, 1.5, 1.5}, {2.0, 2.0, 1.0}, {3.0, 2.5, 1.0}};
  LOG_DEBUG("prediction: particles "
            << 0 << " and " << 1
            << " overlap with dist = " << 0.866025403784439);
  BOOST_CHECK(check_local_dist(pos_0, l, near_min2, near_max2, nnear_min2,
                               nnear_max2) == false);

  const std::vector<Vec3> pos_1 = {
      {2.0, 2.0, 1.0}, {3.0, 2.5, 1.0}, {4.0, 3.0, 1.0}};
  LOG_DEBUG("prediction: particles "
            << 0 << " and " << 2
            << " overlap with dist = " << 2.23606797749979);
  BOOST_CHECK(check_local_dist(pos_1, l, near_min2, near_max2, nnear_min2,
                               nnear_max2) == false);

  const std::vector<Vec3> pos_2 = {
      {0.0, 0.0, 0.0}, {1.17, 0.0, 0.0}, {1.17, 1.0, 0.0}};
  BOOST_CHECK(check_local_dist(pos_2, l, near_min2, near_max2, nnear_min2,
                               nnear_max2) == true);

  const std::vector<Vec3> pos_3 = {
      {1.0, 1.0, 1.0}, {1.0, 2.17, 1.0}, {1.0, 2.17, 2.17}};
  BOOST_CHECK(check_local_dist(pos_3, l, near_min2, near_max2, nnear_min2,
                               nnear_max2) == true);
}

BOOST_AUTO_TEST_CASE(check_snapshot) {
  const std::string snapshot_name = "snapshot_test.h5";

  std::vector<Vec3> pos = {
      {6.0, 10.0, 0.0}, {2.5, 15.0, 3.0},  {9.0, 10.0, 0.0}, {13.0, 11.0, 10.0},
      {7.0, 6.0, 6.5},  {15.0, 8.0, 9.0},  {8.0, 10.0, 0.0}, {9.5, 8.0, 1.3},
      {1.0, 8.0, 4.2},  {17.5, 15.0, 3.0}, {9.5, 10.0, 0.0}};

  std::vector<double> s_bias = {5.5, 3.3, 1.4, 2.5,  16.6, 12.0,
                                9.8, 9.2, 0.8, 22.4, 30.7};

  const unsigned int seed = 40;
  Random mt(seed);

  UpdateConfig update_config;
  update_config.config = 10;

  write_snapshot(snapshot_name, pos, s_bias, mt, update_config);

  std::vector<Vec3> output_pos;
  std::vector<double> output_s_bias;
  Random output_mt(seed + 1);
  UpdateConfig output_update_config;

  read_snapshot(snapshot_name, output_pos, output_s_bias, output_mt,
                output_update_config);

  BOOST_CHECK_EQUAL_COLLECTIONS(output_pos.begin(), output_pos.end(),
                                pos.begin(), pos.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(output_s_bias.begin(), output_s_bias.end(),
                                s_bias.begin(), s_bias.end());
  BOOST_CHECK_EQUAL(mt, output_mt);
  BOOST_CHECK_EQUAL(update_config.config, output_update_config.config);
}

BOOST_AUTO_TEST_CASE(check_dist_between_nonlocal_beads) {
  std::vector<Vec3> pos = {{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0},
                           {1.0, 1.0, 0.0}, {1.0, 2.0, 0.0}, {1.0, 3.0, 0.0},
                           {1.0, 3.0, 1.0}, {1.0, 2.0, 1.0}, {2.0, 2.0, 1.0},
                           {2.0, 2.0, 0.0}, {3.0, 2.0, 0.0}};

  const NonlocalBonds nonlocal_bonds{{{0, 3, 1.5}, {4, 7, 1.5}}};

  const unsigned int nbonds = nonlocal_bonds.get_nbonds();
  std::vector<double> dist(nbonds);

  const double l = 8.0;

  dist_between_nonlocal_beads(pos, l, nonlocal_bonds, dist);

  BOOST_CHECK_CLOSE_FRACTION(dist[0], 1.0, 1e-10);
  BOOST_CHECK_CLOSE_FRACTION(dist[1], 1.0, 1e-10);

  pos.clear();

  pos = {{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
         {1.0, 2.0, 0.0}, {1.0, 3.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, 2.0, 0.0},
         {0.0, 2.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

  dist_between_nonlocal_beads(pos, l, nonlocal_bonds, dist);

  BOOST_CHECK_CLOSE_FRACTION(dist[0], 1.0, 1e-10);
  BOOST_CHECK_CLOSE_FRACTION(dist[1], 1.0, 1e-10);
}
