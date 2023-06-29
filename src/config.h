// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// config.h contains commands related to configurations such as breaking a
// bond, and writes the configuration to the output file

#ifndef HYBRIDMC_CONFIG_H
#define HYBRIDMC_CONFIG_H

#include "json.hpp"
#include <H5Cpp.h>
#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>
#include <vector>

using Config = uint64_t;
using ConfigInt = std::vector<uint64_t>;

// TODO: has to change to account for rc in k index

// TODO: when reading in the pairs, store rc squared rather than rc

class NonlocalBonds {
public:
  using Pairs = std::vector<std::tuple<unsigned int, unsigned int, double>>;

  NonlocalBonds() = default;
  NonlocalBonds(const Pairs &ij);

  // determines if a bond between beads i and j will form by matching the
  // i and j provided to the pairs of beads specified in the json file
  Config get_bond_mask(unsigned int i, unsigned int j) const {
    assert(j > i);

    // check that i and j match with i and j from the json file
    const auto it = std::find(ij_.begin(), ij_.end(), std::make_tuple(i, j));

    if (it == ij_.end())
      return 0;

    // check that the number of configurations do not exceed 64 bits
    assert(ij_.size() <= std::numeric_limits<Config>::digits);

    // convert 1 from 32 to 64 bits and determine position of new bond
    // in bond pattern (using left shift operator)
    return Config(1) << (it - ij_.begin());
  }

  unsigned int get_nbonds() const;

  unsigned int count_bonds() const;

  void write_hdf5(H5::H5Object &obj, const std::string &name) const;

private:
  Pairs ij_;
};

void from_json(const nlohmann::json &json, NonlocalBonds &transient_bonds);
void from_json(const nlohmann::json &json, NonlocalBonds &permanent_bonds);

struct UpdateConfig {
  // integer representing bonding pattern
  Config config;

  // assume that at t = 0, there are no bonds
  UpdateConfig() : config(0) {}

  bool bonded(Config mask) const { return config & mask; }

  bool non_bonded(Config mask) const { return ~config & mask; }

  void flip_bond(Config mask) { config ^= mask; }

  unsigned int count_bonds() const { return std::popcount(config); }
};

struct CountBond {
  unsigned int formed;
  unsigned int broken;
};

class UpdateConfigWriter {
public:
  UpdateConfigWriter(H5::H5File &file);
  void append();

  // vector of the integers of configurations
  ConfigInt config_int;

  void clear() { config_int.clear(); }

private:
  hsize_t file_dims[1];
  H5::DataSpace file_space;
  H5::DataSet dataset_int;
};

#endif
