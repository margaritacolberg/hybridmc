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
#include <iostream>

using Config = uint64_t;
using ConfigInt = std::vector<uint64_t>;

class NonlocalBonds {
public:
  using Pairs = std::vector<std::tuple<unsigned int, unsigned int, double>>;

  NonlocalBonds() = default;
  NonlocalBonds(const Pairs &ij);

  // obtain the number of bonds in the given list of bonds
  unsigned int get_nbonds() const;

  // obtain the bond_index'th bond in the master bonds list
  std::tuple<unsigned int, unsigned int, double> getBond(unsigned int bond_index) const;

  // bond mask is obtained -- a configuration and the rc value packaged in a tuple
  std::tuple<Config, double> get_bond_mask(unsigned int i, unsigned int j) const;

  void write_hdf5(H5::H5Object &obj, const std::string &name) const;

  // print out bonds list
  void printBonds();

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
