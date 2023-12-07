// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// main.cc simulates the dynamics of protein folding using a coarse-grained
// model; each amino acid in the protein is represented by a bead, and the
// beads are connected by local and nonlocal bonds; the dynamics is
// event-driven, and the program converges when each state of the protein has
// been visited roughly equally, giving the final value of the entropy of each
// state;

#include "main_helpers.h"
double max_time = std::numeric_limits<double>::max();

void initialize_pos(System &sys, Random &mt, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    std::optional<std::string> input_name,
                    std::optional<std::string> snapshot_name,
                    const unsigned int t_bonds) {

  if (snapshot_name && std::filesystem::exists(*snapshot_name)) {
    // overwrite existing entries in pos and s_bias vectors with read-in values
    // from hdf5 file
    std::cout << " Reading in snapshot " << *snapshot_name << std::endl;
    read_snapshot(*snapshot_name, sys.pos, sys.s_bias, mt, update_config);
  } else if (input_name && std::filesystem::exists(*input_name)) {
    std::cout << "Reading in input file " << *input_name << std::endl;
    read_input(*input_name, sys.pos);
    init_update_config(sys.pos, update_config, box, p.transient_bonds);
    init_s(sys.s_bias, t_bonds);
  } else {

    // set flag for random initialization success
    bool found = false;
    // set the maximum tries for a random initialization to 10
    int max_init_count = 10;
    // while random initialization not found and 10 attempts not hit, try random init
    while (found == false && max_init_count)
    {
        found = init_pos(sys.pos, box, mt, p);
        max_init_count--;
    };

    // do linear draw if above procedure failed

    if (found == false) {

        std::cout << " Calling linear chain draw." << std::endl;
        draw_linear_chain(sys.pos,p);

    }

    init_s(sys.s_bias, t_bonds);
  }
}

void initialize_system(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, Cells &cells,
                       EventQueue &event_queue) {

  LOG_DEBUG("New step initialization");
  if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                        p.nnear_max2)) {
    throw std::runtime_error("local beads overlap in initialize_system.");
  }

  if (!check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                           p.transient_bonds, p.permanent_bonds)) {
    throw std::runtime_error("nonlocal beads overlap in initialize_system");
  }

  if (cells.ncell < 4) {
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

  if (cells.lcell < check_rc) {
    std::cout << cells.lcell << " is cells.lcell and the check_rc is " << check_rc << std::endl;
    float req_length = cells.ncell * check_rc + 0.0001; // NOLINT(cppcoreguidelines-narrowing-conversions)
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
              Cells &cells, EventQueue &event_queue,
              const unsigned int step, double del_t) {
  LOG_DEBUG("step = " << step);

  double zero_time = 0.0;
  // the current time interval the events are occurring in
  double step_time = step * del_t;

//  std::cout << " At start of step " << step << " local clock time is " << sys.times[0] << std::endl;

  // while events are occurring in step_time,
  while (!event_queue.empty()) {
    // access the minimum time to the next collision, and the indices and
    // collision counters of the beads associated with this collision
    const Event event = event_queue.top();

    if (std::visit(
            [=](auto &&ev) {
              // check for monotonically increasing event times
              assert(ev.t >= zero_time);
              return ev.t > step_time;
            },
            event))
      break;

    event_queue.pop();

    // process collision or cell crossing event
    std::visit(
        [&](auto &&ev) {
          zero_time = ev.t;
          LOG_DEBUG("wall time " << zero_time << " Queue Size is " << event_queue.size());
          process_event(ev, sys, p, box, event_queue, cells, update_config,
                        count_bond);
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

}

void run_trajectory_eq(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, CountBond &count_bond,
                       unsigned int iter,
                       DistWriter &dist_writer, std::vector<double> &dist) {

  LOG_DEBUG("run_trajectory_eq");
  for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
       step++) {
    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

    LOG_DEBUG("max_time = " << max_time);

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t);

    dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
    if (sys.distanceWrite) dist_writer.append(dist);

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
                    CountBond &count_bond,
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

    //set max time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

/*     std::cout << " step = " << step
        << " max_time = " << max_time << " with nsteps = " << p.nsteps
        << std::endl;*/

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    // to check energy conservation
    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t);
    if (sys.distanceWrite)
    {
        dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
        dist_writer.append(dist);
    }

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
    if (E_diff >= 1e-6) {
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
                         CountBond &count_bond,
                         unsigned int iter_wl,
                         bool record_dists,
                         std::vector<double>* dist,
                         DistWriter* dist_writer) {

  LOG_DEBUG("run_trajectory_wl");

  for (unsigned int step = iter_wl * p.nsteps_wl; step < (iter_wl + 1) * p.nsteps_wl; step++) {

    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time
    if (step != 0) {max_time = (step * p.del_t_wl) + 0.001;}

/*    std::cout << " In run_trajectory_wl at iter_wl = " << iter_wl << " step = " << step
        << " max_time = " << max_time << " with nsteps_wl = " << p.nsteps_wl
        << std::endl;*/

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t_wl);

    if (record_dists) {
      dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, *dist);
      dist_writer->append(*dist);
    }

    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));

    if (E_diff >= 1e-6) {
      std::cout << E_diff << " energy difference" << std::endl;
      throw std::runtime_error("energy is not conserved");
    }

  }

  return update_config.config;
}

// Wang-Landau algorithm for estimating entropy
void wang_landau_process(System &sys, Random &mt, const Param &p, const Box &box,
                 UpdateConfig &update_config, CountBond &count_bond,
                 const unsigned int nstates, std::vector<double> &s_bias,
                 DistWriter &dist_writer, std::vector<double> &dist) {

  // amount by which entropy is adjusted
  double gamma = p.gamma;
  // initialize iteration counter at 0
  unsigned int iter_wl = 0;
  unsigned int native_ind = nstates - 1;
  bool record_dists = true;

  while (gamma > p.gamma_f) {
    // iterate over the 2 states
    for (unsigned int i = 0; i < nstates; i++) {
      // run trajectory to get final state
      Config state = run_trajectory_wl(sys, mt, p, box, update_config,
                                       count_bond, iter_wl, record_dists, &dist, &dist_writer);

      if (state == 0) {
        s_bias[native_ind] -= gamma;
      } else {
        s_bias[native_ind] += gamma;
      }

      // update the iteration counter
      iter_wl += 1;
    }
    // normalized gamma updated
    gamma = double(nstates) / double(iter_wl);
  }
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
    p.total_iter = json["total_iter"];
    p.total_iter_eq = json["total_iter_eq"];
    p.pos_scale = json["pos_scale"];
    p.neg_scale = json["neg_scale"];
    p.sig_level = json["sig_level"];
    p.max_nbonds = json["max_nbonds"];
    p.max_g_test_count = json["max_g_test_count"];
    p.flip_req = json["flip_req"];
    p.fail_max = json["fail_max"];
    p.req_dists = json["req_dists"];
}