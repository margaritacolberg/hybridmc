// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// cells.h contains the variables needed to describe a grid of cells that a 3D
// box with periodic boundary conditions is divided into; each cell may contain
// some portion of the protein enclosed in the box

#ifndef HYBRIDMC_CELLS_H
#define HYBRIDMC_CELLS_H

#include <vector>

struct Cells {
  // the inner vector stores the indices of the hard spheres and the outer
  // vector stores the indices of the cells the hard spheres are located in
  std::vector<std::vector<unsigned int>> cells;

  // store the 3D indices of the old cell
  std::vector<unsigned int> cell_x;
  std::vector<unsigned int> cell_y;
  std::vector<unsigned int> cell_z;

  // number of cells in each dimension of the box
  const unsigned int ncell;
  // length, and inverse length of each cell in the box
  const double lcell, inv_lcell;
  Cells(const unsigned int ncell, const double lcell)
      : cells(ncell * ncell * ncell), ncell(ncell), lcell(lcell),
        inv_lcell(1 / lcell) {}
};

#endif
