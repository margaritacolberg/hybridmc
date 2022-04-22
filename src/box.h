// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// box.h implements periodic boundary conditions of a box containing hard
// sphere protein to mimic an infinitely large system

#ifndef HYBRIDMC_BOX_H
#define HYBRIDMC_BOX_H

#include <cmath>

struct Box {
  double l, inv_l;
  Box(double l) : l(l), inv_l(1.0 / l) {}

  void minpos(double &x, double &y, double &z) const {
    x -= std::floor(x * inv_l) * l;
    y -= std::floor(y * inv_l) * l;
    z -= std::floor(z * inv_l) * l;
  }

  void mindist(double &dx, double &dy, double &dz) const {
    dx -= std::round(dx * inv_l) * l;
    dy -= std::round(dy * inv_l) * l;
    dz -= std::round(dz * inv_l) * l;
  }

  // note: x*inv_l represents the x-coordinate in multiples of box lengths
  // (as in, how many times the length of the box was traversed); for
  // example, if x = 7.5 and l = 3, then x is 2.5 boxes; std::floor rounds
  // the number 2.5 to 2, since we want to count full boxes traversed; we
  // then multiply 2 by l to get 6, which is the total distance of 2 boxes
  // traversed completely; so overall, the bead whose x-coordinate at 7.5 is
  // located at x = 7.5 - 6 = 1.5
  //
  // mindist proceeds in the same way as minpos, except here, instead of the
  // result being normalized between 0 and l, it is now between -l/2 and l/2;
  // rounding is done up or down to determine the smallest distance between
  // two beads
};

#endif
