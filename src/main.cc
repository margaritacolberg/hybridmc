// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// main.cc simulates the dynamics of protein folding using a coarse-grained
// model; each amino acid in the protein is represented by a bead, and the
// beads are connected by local and nonlocal bonds; the dynamics is
// event-driven, and the program converges when each state of the protein has
// been visited roughly equally, giving the final value of the entropy of each
// state;

#include "config.h"
#include "crankshaft.h"
#include "entropy.h"
#include "hardspheres.h"
#include "json.hpp"
#include "snapshot.h"
#include "system.h"
#include "writer.h"
#include <H5Cpp.h>
#include <boost/program_options.hpp>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  // restart program from point of interruption
  // see Prus, 2002, Boost.Program_options doc
  //
  // declare the supported options
  po::options_description desc("allowed options");
  desc.add_options()("help,h", "produce help message")(
      "json-file", po::value<std::string>()->required(), "json")(
      "output-file", po::value<std::string>()->required(),
      "hdf5 output")("input-file", po::value<std::string>(), "hdf5 input")(
      "snapshot-file", po::value<std::string>(), "hdf5 snapshot");

  po::positional_options_description pod;
  pod.add("json-file", 1);
  pod.add("output-file", 1);

  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pod).run(),
        vm);

    if (vm.count("help")) {
      std::cout << "usage: hybridmc [options]... json-file output-file"
                << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 0;
    }

    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "hybridmc: " << e.what() << std::endl;
    std::cerr << "try hybridmc --help" << std::endl;
    return 1;
  }

  std::optional<std::string> input_name;
  if (vm.count("input-file")) {
    input_name = vm["input-file"].as<std::string>();
  }

  std::optional<std::string> snapshot_name;
  if (vm.count("snapshot-file")) {
    snapshot_name = vm["snapshot-file"].as<std::string>();
  }

  const std::string json_name = vm["json-file"].as<std::string>();
  const std::string output_name = vm["output-file"].as<std::string>();

  std::cout << "git commit " << VERSION << std::endl;

  std::ifstream input(json_name);
  nlohmann::json json;
  input >> json;

  const Param p = json;
  const Box box{p.length};
  UpdateConfig update_config;

  const unsigned int t_bonds = p.transient_bonds.get_nbonds();
  const unsigned int nbonds = p.nonlocal_bonds.get_nbonds();
  const unsigned int nstates = std::pow(2, t_bonds);
  ConfigInt store_config_int;
  std::vector<uint64_t> config_count(nstates);
  std::vector<double> dist(nbonds);

  // initialize both members of CountBond struct to 0
  CountBond count_bond = {};

  System sys(p.nbeads);

  std::seed_seq seq(p.seeds.begin(), p.seeds.end());
  Random mt(seq);

  if (snapshot_name && std::filesystem::exists(*snapshot_name)) {
    // overwrite existing entries in pos and s_bias vectors with read-in
    // values from hdf5 file
    read_snapshot(*snapshot_name, sys.pos, sys.s_bias, mt, update_config);
  } else if (input_name && std::filesystem::exists(*input_name)) {
    read_input(*input_name, sys.pos);
    init_update_config(sys.pos, update_config, box, p.rc2, p.transient_bonds);
    init_s(sys.s_bias, t_bonds);
  } else {
    init_pos(sys.pos, box, mt, p);
    init_s(sys.s_bias, t_bonds);
  }

  if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                        p.nnear_max2)) {
    throw std::runtime_error("local beads overlap");
  }

  if (!check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.stair2, p.p_rc2,
                           p.transient_bonds, p.permanent_bonds)) {
    throw std::runtime_error("nonlocal beads overlap");
  }

  // open output file
  H5::H5File file(output_name + ".tmp", H5F_ACC_TRUNC);

  UpdateConfigWriter update_config_writer(file);

  H5::Group group{file.createGroup("unique_config")};
  PosWriter pos_writer{group, "pos", p.nbeads};
  VelWriter vel_writer{group, "vel", p.nbeads};
  ConfigWriter config_writer{group, "config"};
  ConfigCountWriter config_count_writer{file, "config_count", nstates};
  DistWriter dist_writer{file, "dist", nbonds};

  config_count_writer.append(config_count);
  std::set<Config> store_config;

  const hsize_t mem_dims[1] = {sys.s_bias.size()};
  H5::DataSpace mem_space(1, mem_dims);
  H5::DataSet dataset_s_bias =
      file.createDataSet("s_bias", H5::PredType::NATIVE_FLOAT, mem_space);

  H5::Attribute attr_sigma =
      file.createAttribute("sigma", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_sigma.write(H5::PredType::NATIVE_DOUBLE, &p.sigma);
  H5::Attribute attr_length =
      file.createAttribute("length", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_length.write(H5::PredType::NATIVE_DOUBLE, &box.l);
  H5::Attribute attr_nsteps =
      file.createAttribute("nsteps", H5::PredType::NATIVE_UINT, H5S_SCALAR);
  attr_nsteps.write(H5::PredType::NATIVE_UINT, &p.nsteps);

  p.nonlocal_bonds.write_hdf5(file, "nonlocal_bonds");
  p.transient_bonds.write_hdf5(file, "transient_bonds");
  p.permanent_bonds.write_hdf5(file, "permanent_bonds");

  unsigned int g_test_count = 0;

  // the BIIIIIIIG loop
  do {
    // reset bead clocks, counters and wall time
    for (unsigned int i = 0; i < p.nbeads; i++) {
      sys.times[i] = 0.0;
      sys.counter[i] = 0.0;
    }

    double wall_time = 0.0;

    // reset configuration count
    store_config_int.clear();
    std::fill(config_count.begin(), config_count.end(), 0);

    for (unsigned int iter = 0; iter < p.total_iter; iter++) {
      Cells cells{p.ncell, p.length / p.ncell};

      if (!(cells.ncell >= 4)) {
        throw std::invalid_argument("bead ncell must be at least 4");
      }

      double check_rc = p.rc;
      if (p.stair) {
        check_rc = *p.stair;
      }

      if (!(cells.lcell >= check_rc)) {
        throw std::invalid_argument("bead lcell must be at least rc");
      }

      // split the box into smaller cells, and store the beads in each
      // cell
      init_cells(sys.pos, box, cells);

      init_vel(sys.vel, mt, p.temp, p.m);

      // to check energy conservation
      const double tot_E_before =
          compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

      EventQueue event_queue;

      // fill priority queue with cell crossings of all particles
      init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                       event_queue, sys.times, cells);

      // fill priority queue with nearest bond events between all
      // particle pairs
      init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                               event_queue, sys.times, p.near_min2,
                               p.near_max2);

      // fill priority queue with next-nearest bond events between all
      // particle pairs
      init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                                event_queue, sys.times, p.nnear_min2,
                                p.nnear_max2);

      // fill priority queue with collisions of all particle pairs
      // (executed once, with an initial prediction that fills the entire
      // priority queue)
      add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.rc2,
                               p.stair2, p.p_rc2, box, sys.counter, event_queue,
                               sys.times, cells, p.transient_bonds,
                               p.permanent_bonds, update_config, p.max_nbonds);

      // assume that the entire time during which the beads are
      // undergoing events can be divided into intervals of length
      // p.del_t; the total number of such intervals is p.nsteps (thus,
      // the variable called step marks the intervals until the
      // p.nsteps-th interval is reached); each iteration of the
      // BIIIIIIIG loop will run until update_time and then dump the
      // output to the hdf5 file
      for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
           step++) {
        LOG_DEBUG("step = " << step);

        // the current time interval the events are occurring in
        double step_time = step * p.del_t;

        // while events are occurring in step_time,
        while (!event_queue.empty()) {
          // access the minimum time to the next collision, and the
          // indices and collision counters of the beads associated
          // with this collision
          const Event event = event_queue.top();

          if (std::visit(
                  [=](auto &&ev) {
                    // check for monotonically increasing event times
                    assert(ev.t >= wall_time);
                    return ev.t > step_time;
                  },
                  event))
            break;

          event_queue.pop();

          // process collision or cell crossing event
          std::visit(
              [&](auto &&ev) {
                wall_time = ev.t;
                LOG_DEBUG("wall time " << wall_time);
                process_event(ev, sys, p, box, event_queue, cells,
                              update_config, count_bond);
              },
              event);
        }

        // update positions at the moment the last collision in p.del_t
        // occurred
        for (unsigned int i = 0; i < p.nbeads; i++) {
          update_pos(sys.pos[i], sys.vel[i], sys.times[i], step_time);
          assert(check_overlap(i, sys.pos, sys.vel, sys.times, p.rh2, box));
        }

        // check bond distances of nearest and next-nearest beads
        assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2,
                                p.nnear_min2, p.nnear_max2));
        // check bond distances of nonlocal beads
        assert(check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.stair2,
                                   p.p_rc2, p.transient_bonds,
                                   p.permanent_bonds));

        dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
        dist_writer.append(dist);

        if (step % p.write_step == 0) {
          // store the integer of the configuration and the time of
          // the event
          store_config_int.emplace_back(update_config.config);
          update_config_writer.config_int.emplace_back(update_config.config);

          // if a transient bond forms, check if configuration has
          // been previously visited by comparing configuration to
          // set of saved configurations; if not visited, save
          // configuration to the set and write positions of beads to
          // file
          if (update_config.config == 1 &&
              store_config.insert(update_config.config).second) {
            pos_writer.append(sys.pos);
            vel_writer.append(sys.vel);
            config_writer.append(update_config.config);
          }
        }

        // update time
        wall_time = step_time;

        const double tot_E_during =
            compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);
        const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));
        if (!(E_diff < 1e-6)) {
          std::cout << E_diff << " energy difference" << std::endl;
          throw std::runtime_error("energy is not conserved");
        }
      }

      for (unsigned int i = 0; i < p.mc_moves; i++) {
        crankshaft(sys.pos, update_config, box, p.near_min2, p.near_max2,
                   p.nnear_min2, p.nnear_max2, p.rc2, p.rh2, p.stair2, p.p_rc2,
                   p.transient_bonds, p.permanent_bonds, mt, sys.s_bias);
        for (unsigned int j = 1; j < p.mc_write; j++) {
          if (i == ((j * p.mc_moves) / p.mc_write)) {
            update_config_writer.config_int.emplace_back(update_config.config);
          }
        }
      }

      assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2,
                              p.nnear_min2, p.nnear_max2));
      assert(check_nonlocal_dist(sys.pos, box, p.rc2, p.rh2, p.stair2, p.p_rc2,
                                 p.transient_bonds, p.permanent_bonds));

      // store the configurations in the hdf5 file
      update_config_writer.append();
      update_config_writer.clear();
    }

    // count the number of times each configuration is visited
    for (unsigned int i = 0; i < store_config_int.size(); i++) {
      config_count[store_config_int[i]]++;
    }

    compute_entropy(config_count, sys.s_bias, nstates, p.pos_scale,
                    p.neg_scale);
    dataset_s_bias.write(&sys.s_bias[0], H5::PredType::NATIVE_DOUBLE,
                         mem_space);

    if (snapshot_name) {
      // write new snapshot to temporary file
      const std::string snapshot_tmp = *snapshot_name + ".tmp";
      write_snapshot(snapshot_tmp, sys.pos, sys.s_bias, mt, update_config);
      // atomically overwrite old snapshot with new snapshot
      std::filesystem::rename(snapshot_tmp, *snapshot_name);
    }

    config_count_writer.append(config_count);

  } while (++g_test_count < p.max_g_test_count &&
           !g_test(config_count, nstates, p.sig_level));

  file.close();

  std::filesystem::rename(output_name + ".tmp", output_name);

  return 0;
}

void from_json(const nlohmann::json &json, Param &p) {
  p.m = json["m"];
  p.sigma = json["sigma_bb"];
  p.sigma2 = p.sigma * p.sigma;
  p.near_min = json["near_min"];
  p.near_max = json["near_max"];
  p.near_min2 = p.near_min * p.near_min;
  p.near_max2 = p.near_max * p.near_max;
  p.nnear_min = json["nnear_min"];
  p.nnear_max = json["nnear_max"];
  p.nnear_min2 = p.nnear_min * p.nnear_min;
  p.nnear_max2 = p.nnear_max * p.nnear_max;
  p.rh = json["rh"];
  p.rc = json["rc"];
  p.rh2 = p.rh * p.rh;
  p.rc2 = p.rc * p.rc;

  if (json.count("stair") != 0) {
    p.stair = json["stair"];
    p.stair2 = *p.stair * *p.stair;
  }

  p.nonlocal_bonds = json["nonlocal_bonds"];
  p.transient_bonds = json["transient_bonds"];
  p.permanent_bonds = json["permanent_bonds"];

  if (json.count("stair_bonds") != 0) {
    p.stair_bonds = json["stair_bonds"];
  }

  if (json.count("p_rc") != 0) {
    p.p_rc = json["p_rc"];
    p.p_rc2 = *p.p_rc * *p.p_rc;
  }

  std::cout << *p.stair << std::endl;
  if (p.stair && ((*p.stair < p.rc) || (*p.stair < *p.p_rc))) {
    throw std::runtime_error("stair boundary is smaller than rc");
  }

  p.tries = json["tries"];
  p.nbeads = json["nbeads"];
  p.length = json["length"];
  p.ncell = json["ncell"];
  p.nsteps = json["nsteps"];
  p.del_t = json["del_t"];
  p.write_step = json["write_step"];
  p.seeds = json["seeds"].get<std::vector<unsigned int>>();
  p.temp = json["temp"];
  p.mc_moves = json["mc_moves"];
  p.mc_write = json["mc_write"];
  p.total_iter = json["total_iter"];
  p.pos_scale = json["pos_scale"];
  p.neg_scale = json["neg_scale"];
  p.sig_level = json["sig_level"];
  p.max_nbonds = json["max_nbonds"];
  p.max_g_test_count = json["max_g_test_count"];
}
