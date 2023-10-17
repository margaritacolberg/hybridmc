// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// main.cc simulates the dynamics of protein folding using a coarse-grained
// model; each amino acid in the protein is represented by a bead, and the
// beads are connected by local and nonlocal bonds; the dynamics is
// event-driven, and the program converges when each state of the protein has
// been visited roughly equally, giving the final value of the entropy of each
// state;

#include "crankshaft.h"
#include "entropy.h"
#include "hardspheres.h"
#include "json.hpp"
#include "snapshot.h"
#include "writer.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <H5Cpp.h>
#include <boost/program_options.hpp>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>

double max_time = std::numeric_limits<double>::max();
namespace po = boost::program_options;
namespace py = pybind11;

void initialize_pos(System &sys, Random &mt, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    std::optional<std::string> input_name,
                    std::optional<std::string> snapshot_name,
                    const unsigned int t_bonds) {
  if (snapshot_name && std::filesystem::exists(*snapshot_name)) {
    // overwrite existing entries in pos and s_bias vectors with read-in values
    // from hdf5 file
    read_snapshot(*snapshot_name, sys.pos, sys.s_bias, mt, update_config);
  } else if (input_name && std::filesystem::exists(*input_name)) {
    read_input(*input_name, sys.pos);
    init_update_config(sys.pos, update_config, box, p.transient_bonds);
    init_s(sys.s_bias, t_bonds);
  } else {
    init_pos(sys.pos, box, mt, p);
    init_s(sys.s_bias, t_bonds);
  }
}

void initialize_system(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, Cells &cells,
                       EventQueue &event_queue) {


  LOG_DEBUG("New step initialization");
  if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                        p.nnear_max2)) {
    throw std::runtime_error("local beads overlap");
  }

  if (!check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                           p.transient_bonds, p.permanent_bonds)) {
    throw std::runtime_error("nonlocal beads overlap");
  }

  if (!(cells.ncell >= 4)) {
    throw std::invalid_argument("bead ncell must be at least 4");
  }

  double check_rc = 0;
  for (unsigned int i=0; i < p.nonlocal_bonds.get_nbonds(); i++) {

      if (p.nonlocal_bonds.getrc(i) > check_rc) {
          check_rc = p.nonlocal_bonds.getrc(i);
      }
  }

  if (p.stair) {
    check_rc = *p.stair;
  }

  if (!(cells.lcell >= check_rc)) {
    std::cout << cells.lcell << " is cells.lcell and the check_rc is " << check_rc << std::endl;
    float req_length = cells.ncell * check_rc + 0.0001;
    std::cout << "Length has to be at least " << req_length << std::endl;
    throw std::invalid_argument("bead lcell must be at least rc. Increase length to at least above value");
  }

  // split the box into smaller cells, and store the beads in each cell
  init_cells(sys.pos, box, cells);

  init_vel(sys.vel, mt, p.temp, p.m);

  // fill priority queue with cell crossings of all particles
  init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter, event_queue,
                   sys.times, cells);

  // fill priority queue with nearest bond events between all particle pairs
  init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                           event_queue, sys.times, p.near_min2, p.near_max2);

  // fill priority queue with next-nearest bond events between all particle
  // pairs
  init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                            event_queue, sys.times, p.nnear_min2, p.nnear_max2);

  // fill priority queue with collisions of all particle pairs (executed once,
  // with an initial prediction that fills the entire priority queue)
  add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.stair2, box, sys.counter, event_queue, sys.times,
                           cells, p.transient_bonds, p.permanent_bonds,
                           update_config, p.max_nbonds);
}

void run_step(System &sys, const Param &p, const Box &box,
              UpdateConfig &update_config, CountBond &count_bond,
              double wall_time, Cells &cells, EventQueue &event_queue,
              const unsigned int step, double del_t) {
  LOG_DEBUG("step = " << step);

  // the current time interval the events are occurring in
  double step_time = step * del_t;

  // while events are occurring in step_time,
  while (!event_queue.empty()) {
    // access the minimum time to the next collision, and the indices and
    // collision counters of the beads associated with this collision
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
          LOG_DEBUG("wall time " << wall_time << " Queue Size is " << event_queue.size());
          process_event(ev, sys, p, box, event_queue, cells, update_config,
                        count_bond);
          //if (wall_time > 15.1) {exit(0);}
        },
        event);
  }

  // update positions at the moment the last collision in del_t occurred
  for (unsigned int i = 0; i < p.nbeads; i++) {
    update_pos(sys.pos[i], sys.vel[i], sys.times[i], step_time);
    assert(check_overlap(i, sys.pos, sys.vel, sys.times, p.rh2, box));
  }

  // check bond distances of nearest and next-nearest beads
  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));
  // check bond distances of nonlocal beads
  assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds));

  // update time
  wall_time = step_time;
}

void run_trajectory_eq(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, CountBond &count_bond,
                       double wall_time, unsigned int iter) {

  LOG_DEBUG("run_trajectory_eq");
  for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
       step++) {
    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time based on wall_time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

    LOG_DEBUG("max_time = " << max_time << " wall_time = " << wall_time);

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    run_step(sys, p, box, update_config, count_bond, wall_time, cells,
             event_queue, step, p.del_t);
  }

  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));
  assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds));
}

void run_trajectory(System &sys, Random &mt, const Param &p, const Box &box,
                    std::vector<double> &dist, UpdateConfig &update_config,
                    UpdateConfigWriter &update_config_writer,
                    PosWriter &pos_writer, VelWriter &vel_writer,
                    ConfigWriter &config_writer, DistWriter &dist_writer,
                    std::set<Config> &store_config, ConfigInt &store_config_int,
                    CountBond &count_bond, double wall_time,
                    unsigned int iter) {

  LOG_DEBUG("run_trajectory");

  // Do Molecular Dynamics moves

  // assume that the entire time during which the beads are undergoing events
  // can be divided into intervals of length p.del_t; the total number of such
  // intervals is p.nsteps (thus, the variable called step marks the intervals
  // until the p.nsteps-th interval is reached); each iteration of the
  // BIIIIIIIG loop will run until update_time and then dump the output to the
  // hdf5 file
  for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
       step++) {
    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time based on wall_time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}
    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    // to check energy conservation
    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, wall_time, cells,
             event_queue, step, p.del_t);

    dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
    dist_writer.append(dist);

    if (step % p.write_step == 0) {
      // store the integer of the configuration and the time of the event

      //store_config_int.emplace_back(update_config.config);
      //update_config_writer.config_int.emplace_back(update_config.config);

      // if a transient bond forms, check if configuration has been previously
      // visited by comparing configuration to set of saved configurations; if
      // not visited, save configuration to the set and write positions of
      // beads to file
      if (update_config.config == 1 &&
          store_config.insert(update_config.config).second) {
        pos_writer.append(sys.pos);
        vel_writer.append(sys.vel);
        config_writer.append(update_config.config);
      }
    }

    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);
    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));
    if (!(E_diff < 1e-6)) {
      std::cout << E_diff << " energy difference" << std::endl;
      throw std::runtime_error("energy is not conserved");
    }
  }

  // Do Monte Carlo moves (mc moves)
  for (unsigned int i = 0; i < p.mc_moves; i++) {
      crankshaft(sys.pos, update_config, box, p.near_min2, p.near_max2,
                 p.nnear_min2, p.nnear_max2, p.rh2, p.stair2,
                 p.transient_bonds, p.permanent_bonds, mt, sys.s_bias);
      /* for (unsigned int j = 1; j < p.mc_write; j++) {
          if (i == ((j * p.mc_moves) / p.mc_write)) {
              update_config_writer.config_int.emplace_back(update_config.config);
          }
      }*/
  }

  // store the integer of the configuration and the time of the event
  store_config_int.emplace_back(update_config.config);
  update_config_writer.config_int.emplace_back(update_config.config);

  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));
  assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds));

  // store the configurations in the hdf5 file
  update_config_writer.append();
  update_config_writer.clear();
}

Config run_trajectory_wl(System &sys, Random &mt, const Param &p,
                         const Box &box, UpdateConfig &update_config,
                         CountBond &count_bond, double wall_time,
                         const unsigned int iter_wl) {

  LOG_DEBUG("run_trajectory_wl");

  for (unsigned int step = iter_wl * p.nsteps_wl;
       step < (iter_wl + 1) * p.nsteps_wl; step++) {

    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time based on wall_time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, wall_time, cells,
             event_queue, step, p.del_t_wl);

    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);
    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));
    if (!(E_diff < 1e-6)) {
      std::cout << E_diff << " energy difference" << std::endl;
      throw std::runtime_error("energy is not conserved");
    }
  }

  return update_config.config;
}

// Wang-Landau algorithm for estimating entropy
void wang_landau(System &sys, Random &mt, const Param &p, const Box &box,
                 UpdateConfig &update_config, CountBond &count_bond,
                 const unsigned int nstates, std::vector<double> &s_bias) {
  // amount by which entropy is adjusted
  double gamma = p.gamma;
  unsigned int iter_wl = 0;
  unsigned int native_ind = nstates - 1;
  double wall_time = 0.0;

  while (gamma > p.gamma_f) {
    // iterate over the 2 states
    for (unsigned int i = 0; i < nstates; i++) {
      // run trajectory to get final state
      Config state = run_trajectory_wl(sys, mt, p, box, update_config,
                                       count_bond, wall_time, iter_wl);

      if (state == 0) {
        s_bias[native_ind] -= gamma;
      } else {
        s_bias[native_ind] += gamma;
      }
    }

        iter_wl += 1;
        gamma = 1.0 / double(iter_wl);
    }
}

// wang landau plugin for python function
std::vector<double> wang_landau_process(std::string json_name, std::optional<std::string> input_name) {

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

    if (input_name.has_value() && std::filesystem::exists(*input_name)) {
        read_input(*input_name, sys.pos);
        init_update_config(sys.pos, update_config, box, p.transient_bonds);
        init_s(sys.s_bias, t_bonds);
    } else {
        init_pos(sys.pos, box, mt, p);
        init_s(sys.s_bias, t_bonds);
    }

    if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2)) {
        throw std::runtime_error("local beads overlap");
    }

    if (!check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds)) {
        throw std::runtime_error("nonlocal beads overlap");
    }

    double gamma = p.gamma;
    //get bead index for transient pair and rc value for it
    std::tuple<unsigned int, unsigned int, double> transient_pair = p.transient_bonds.getBond(0);
    int bead1 = std::get<0>(transient_pair);
    int bead2 = std::get<1>(transient_pair);
   // std::cout << "bead1: " << bead1 << " bead2: " << bead2 << std::endl;
   // the minimum rc2 can be is the target bonding rc for this nonlocal bond; get this value and assign to rc_min2
    double rc_min2 = std::get<2>(transient_pair);
    // take square root to get the minimum bonding distance value
    double rc_min = std::sqrt(rc_min2);
    unsigned int iter_wl = 0;
    unsigned int native_ind = nstates - 1;
    double wall_time = 0.0;
    std::vector<double> distance_values;

    while (gamma > p.gamma_f) {
        // iterate over the 2 states
        for (unsigned int i = 0; i < nstates; i++) {
            // run trajectory to get final state
            Config state = run_trajectory_wl(sys, mt, p, box, update_config,
                                             count_bond, wall_time, iter_wl);

            if (state == 0) {
                sys.s_bias[native_ind] -= gamma;
            } else {
                sys.s_bias[native_ind] += gamma;
            }

            double dx = sys.pos[bead2].x - sys.pos[bead1].x;
            double dy = sys.pos[bead2].y - sys.pos[bead1].y;
            double dz = sys.pos[bead2].z - sys.pos[bead1].z;
            box.mindist(dx, dy, dz);

            const double dist = sqrt(dx * dx + dy * dy + dz * dz);

            //std::cout << " recording everything" << dist << " sb = " << sys.s_bias[native_ind] << std::endl;

            if (dist > rc_min){
                //std::cout << " recording " << dist << " sb = " << sys.s_bias[native_ind] << std::endl;
                distance_values.push_back(dist);
            }
        }

        iter_wl += 1;
        gamma = 1.0 / double(iter_wl);
    }

    //sort distance vector
    std::sort(distance_values.begin(), distance_values.end());

    std::cout << "sbias is " << sys.s_bias[0] - sys.s_bias[1] << std::endl;
    std::cout << "native sbias is " << sys.s_bias[1] << std::endl;
    std::cout << "non native sbias is " << sys.s_bias[0] << std::endl;

    // create return vector
    std::vector<double> return_info;
    // get s bias value
    return_info.push_back(sys.s_bias[0] - sys.s_bias[1]);
    // get rc_min
    return_info.push_back(rc_min);
    // get size of distance value array
    int n_rcs = distance_values.size();
    // want integer index of n_rcs of 0.10 (10%), 0.30 (30%) etc
    double rc_middle = distance_values[int(n_rcs * 0.005)];
    double rc_outermost = distance_values[int(n_rcs * 0.05 )];

    std::cout << "n_rcs: " << n_rcs << " rc_middle: " << rc_middle << " rc_outermost: " << rc_outermost << std::endl;
    std::cout << "True rc " << rc_min << " rc_min: " << distance_values[0] << " rc_max: " << distance_values[n_rcs - 1] << std::endl;

    return_info.push_back(rc_middle);
    return_info.push_back(rc_outermost);

  return return_info;
}

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
  Param p_eq = p;
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

  // equilibrate initial structure
  //
  // reset bead clocks, counters, and wall time
  for (unsigned int i = 0; i < p.nbeads; i++) {
    sys.times[i] = 0.0;
    sys.counter[i] = 0.0;
  }

  double wall_time = 0.0;
  UpdateConfig update_config_eq;

  p_eq.nsteps = p.nsteps_eq;
  p_eq.total_iter = p.total_iter_eq;

  initialize_pos(sys, mt, p_eq, box, update_config_eq, input_name,
                 snapshot_name, t_bonds);

  for (unsigned int iter = 0; iter < p_eq.total_iter; iter++) {
    run_trajectory_eq(sys, mt, p_eq, box, update_config_eq, count_bond,
                      wall_time, iter);
  }

  // Wang-Landau
  init_update_config(sys.pos, update_config, box, p.transient_bonds);

  for (unsigned int i = 0; i < p.nbeads; i++) {
    sys.times[i] = 0.0;
    sys.counter[i] = 0.0;
  }

  wang_landau(sys, mt, p, box, update_config, count_bond, nstates, sys.s_bias);

  const hsize_t mem_dims[1] = {sys.s_bias.size()};
  H5::DataSpace mem_space(1, mem_dims);
  H5::DataSet dataset_s_bias =
      file.createDataSet("s_bias", H5::PredType::NATIVE_FLOAT, mem_space);

  // parameters for checking if g test is validated or not
  unsigned int g_test_count = 0;
  double flipping_rate = 0;
  bool done_flip = false;
  bool done_g_test = false;

  // the BIIIIIIIG loop
  do {
    // reset bead clocks, counters and wall time
    for (unsigned int i = 0; i < p.nbeads; i++) {
      sys.times[i] = 0.0;
      sys.counter[i] = 0.0;
    }

    wall_time = 0.0;

    // reset configuration count
    store_config_int.clear();
    std::fill(config_count.begin(), config_count.end(), 0);

    // add this to see to only count bond events in each iteration
    count_bond.formed = 0;
    count_bond.broken = 0;

    for (unsigned int iter = 0; iter < p.total_iter; iter++) {
      run_trajectory(sys, mt, p, box, dist, update_config, update_config_writer,
                     pos_writer, vel_writer, config_writer, dist_writer,
                     store_config, store_config_int, count_bond, wall_time,
                     iter);
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

    double flips = double(count_bond.formed + count_bond.broken);
    auto  stateCount = std::reduce(config_count.begin(), config_count.end());
    flipping_rate = flips / stateCount;

    done_g_test = g_test(config_count,nstates, p.sig_level);
    if (flipping_rate > p.flip_req) done_flip = true;

    std::cout << "Working on output_file: " << output_name << std::endl;
    std::cout << " In iteration " << g_test_count << " stateCount = " << stateCount << " flips = " << flips << " flip rate = " << flipping_rate
              << " must be greater than " << p.flip_req << " done_flip = " << done_flip << " done_g = " << done_g_test << std::endl;
    g_test_count++;
    if (g_test_count >= p.max_g_test_count)
    {
	    std::cout << "  Done the maximum number of iterations: " << g_test_count <<  " without convergence." << std::endl;
	    break;
    }
  } while (!done_g_test or !done_flip);

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
    p.rh2 = p.rh * p.rh;

    // if stair var is not false aka 0
    if (json.count("stair") != 0) {
        p.stair = json["stair"];
        p.stair2 = *p.stair * *p.stair;
    }

    p.nonlocal_bonds = json["nonlocal_bonds"];
    p.transient_bonds = json["transient_bonds"];
    p.permanent_bonds = json["permanent_bonds"];

    // define the stair bonds
    /*if (json.count("stair_bonds") != 0) {
        p.stair_bonds = json["stair_bonds"];
    }*/

    p.tries = json["tries"];
    p.nbeads = json["nbeads"];
    p.length = json["length"];
    p.ncell = json["ncell"];
    p.nsteps = json["nsteps"];
    p.nsteps_eq = json["nsteps_eq"];
    p.del_t = json["del_t"];
    p.nsteps_wl = json["nsteps_wl"];
    p.del_t_wl = json["del_t_wl"];
    p.gamma = json["gamma"];
    p.gamma_f = json["gamma_f"];
    p.write_step = json["write_step"];
    p.seeds = json["seeds"].get<std::vector<unsigned int>>();
    p.temp = json["temp"];
    p.mc_moves = json["mc_moves"];
    p.mc_write = json["mc_write"];
    p.total_iter = json["total_iter"];
    p.total_iter_eq = json["total_iter_eq"];
    p.pos_scale = json["pos_scale"];
    p.neg_scale = json["neg_scale"];
    p.sig_level = json["sig_level"];
    p.max_nbonds = json["max_nbonds"];
    p.max_g_test_count = json["max_g_test_count"];
    p.flip_req = json["flip_req"];
}

// create pybind11 module for wang_landau function
PYBIND11_MODULE(wang_landau, m) {
    m.doc() = "pybind11 hybridmc plugin for wang_landau function";

    using namespace pybind11::literals;
    m.def("WL_process", &wang_landau_process, "A function that conducts the wang landau algorithm",
          "json_name"_a, "h5_input_name"_a=py::none());

}
