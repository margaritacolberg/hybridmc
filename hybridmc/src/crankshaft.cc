// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// crankshaft.cc performs crankshaft Monte Carlo rotations along randomly
// chosen bonds to create a new configuration of the protein

#include "crankshaft.h"
#include <algorithm>
#include <cmath>

//#define LOCAL_DEBUG
//#define LOCAL_VERBOSE

void swapMC(System &sys, Random &mt, const Box &box, const Param &p) {
    int numEnsemble = sys.ensemble.size();
    if (numEnsemble < 1) return; // no move possible

    std::vector<Vec3> pos = sys.pos;
    Molecule swapMolecule( pos.size() );
    swapMolecule.setPositions(pos); // converts Vec3 positions into vector<double> array
    UpdateConfig current_config = config_int(pos, box, p.transient_bonds);
    int current_state = current_config.config;

#ifdef LOCAL_VERY_VERBOSE
    std::cout << " Before swap, state is: "  << current_state << std::endl;
    swapMolecule.printPositions(-1);
    std::cout << " Ensemble is:" << std::endl;
    printEnsemble(sys.ensemble);
#endif

    std::uniform_int_distribution<> index_probability(0, numEnsemble-1);

    int s_index = index_probability(mt);
    std::vector<Vec3> trial_pos = sys.ensemble[s_index].getVec3Positions();
    UpdateConfig trial_config = config_int(trial_pos, box, p.transient_bonds);
    int trial_state = trial_config.config;


#ifdef LOCAL_VERY_VERBOSE
    std::cout << "  Trying to swap with ensemble index " << s_index << " which has state "
            << trial_state << std::endl;
#endif
    //  replace ensemble state if not the same state
    if (trial_state != current_state )
    {
        sys.ensemble[s_index]  = swapMolecule; // replace incorrect state in ensemble with current state
#ifdef LOCAL_VERBOSE
        std::cout << "  Since state was different, swapping out ensemble index " << s_index << " which had state "
            << trial_state << " replaced with current state and return." << std::endl;
#endif
        return;
    }

    std::swap(sys.ensemble[s_index], swapMolecule);

    //std::cout << std::endl << " After swap of molecule " << s_index << " swap molecule is:" << std::endl;
    //swapMolecule.printPositions(-1);
    //std::cout << " Ensemble is:" << std::endl;
    //printEnsemble(sys.ensemble);

    sys.pos = swapMolecule.getVec3Positions();

}

void rodrigues_rotation(const std::vector<Vec3> &pos, const double theta,
                        std::vector<Vec3> &pos_trial, const unsigned int ind,
                        const Box &box) {
  const unsigned int nbeads = pos.size();
  assert(pos_trial.size() == nbeads);

  assert(1 <= ind && ind <= nbeads - 2);
  assert(0 <= theta && theta <= 2 * M_PI);

  LOG_DEBUG("Rodrigues rotation of particle " << ind << " along angle "
                                              << theta);

  // difference between the positions of the chosen bead in [1, nbeads - 1],
  // and its neighbor in [0, nbeads - 2]
  double dx = pos[ind].x - pos[ind - 1].x;
  double dy = pos[ind].y - pos[ind - 1].y;
  double dz = pos[ind].z - pos[ind - 1].z;
  box.mindist(dx, dy, dz);

  double inv_diff_norm = 1 / std::sqrt(dx * dx + dy * dy + dz * dz);

  // unit vector defining a rotation axis (the rotation axis is along the
  // bond between ind and ind - 1)
  double k_x = dx * inv_diff_norm;
  double k_y = dy * inv_diff_norm;
  double k_z = dz * inv_diff_norm;

  const double sin_theta = std::sin(theta);
  const double cos_theta = 1 - std::cos(theta);

  // Rodrigues rotation matrix
  // see: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
  const double Rxx = 1 + cos_theta * (-k_z * k_z - k_y * k_y);
  const double Rxy = -sin_theta * k_z + cos_theta * k_y * k_x;
  const double Rxz = sin_theta * k_y + cos_theta * k_z * k_x;
  const double Ryx = sin_theta * k_z + cos_theta * k_x * k_y;
  const double Ryy = 1 + cos_theta * (-k_z * k_z - k_x * k_x);
  const double Ryz = -sin_theta * k_x + cos_theta * k_z * k_y;
  const double Rzx = -sin_theta * k_y + cos_theta * k_x * k_z;
  const double Rzy = sin_theta * k_x + cos_theta * k_y * k_z;
  const double Rzz = 1 + cos_theta * (-k_y * k_y - k_x * k_x);

  std::copy(pos.begin(), pos.begin() + ind + 1, pos_trial.begin());

  for (unsigned int i = ind + 1; i < nbeads; i++) {
    double old_diff_x = pos[i].x - pos[ind].x;
    double old_diff_y = pos[i].y - pos[ind].y;
    double old_diff_z = pos[i].z - pos[ind].z;
    box.mindist(old_diff_x, old_diff_y, old_diff_z);

    // rotate the bond between ind and ind - 1
    double new_diff_x = Rxx * old_diff_x + Rxy * old_diff_y + Rxz * old_diff_z;
    double new_diff_y = Ryx * old_diff_x + Ryy * old_diff_y + Ryz * old_diff_z;
    double new_diff_z = Rzx * old_diff_x + Rzy * old_diff_y + Rzz * old_diff_z;
    box.mindist(new_diff_x, new_diff_y, new_diff_z);

    // update the positions of the beads affected by the rotation (position
    // of beads 0 to ind remain unchanged)
    pos_trial[i].x = pos[ind].x + new_diff_x;
    pos_trial[i].y = pos[ind].y + new_diff_y;
    pos_trial[i].z = pos[ind].z + new_diff_z;
  }
}

bool check_bond_if_crankshaft(const double min_dist2, const double max_dist2,
                              double dx, double dy, double dz) {
  double dist2 = dx * dx + dy * dy + dz * dz;

  // ! accounts for comparison with NaN values
  if (!(dist2 > min_dist2 && dist2 < max_dist2)) {
    LOG_DEBUG("check bond distance " << std::sqrt(dist2));
    return false;
  }

  return true;
}

bool check_local_dist_if_crankshaft(const std::vector<Vec3> &pos_trial,
                                    const Box &box, const double near_min2,
                                    const double near_max2,
                                    const double nnear_min2,
                                    const double nnear_max2) {
  const unsigned int nbeads = pos_trial.size();

  // nearest neighbor check
  for (unsigned int i = 0; i + 1 < nbeads; i++) {
    double dx = pos_trial[i + 1].x - pos_trial[i].x;
    double dy = pos_trial[i + 1].y - pos_trial[i].y;
    double dz = pos_trial[i + 1].z - pos_trial[i].z;
    box.mindist(dx, dy, dz);

    if (!check_bond_if_crankshaft(near_min2, near_max2, dx, dy, dz)) {
      LOG_DEBUG("nearest neighbors " << i << " and " << i + 1 << " overlap");
      return false;
    }
  }

  // next-nearest neighbor check
  for (unsigned int i = 0; i + 2 < nbeads; i++) {
    double dx = pos_trial[i + 2].x - pos_trial[i].x;
    double dy = pos_trial[i + 2].y - pos_trial[i].y;
    double dz = pos_trial[i + 2].z - pos_trial[i].z;
    box.mindist(dx, dy, dz);

    if (!check_bond_if_crankshaft(nnear_min2, nnear_max2, dx, dy, dz)) {
      LOG_DEBUG("next-nearest neighbors " << i << " and " << i + 2
                                          << " overlap");
      return false;
    }
  }

  return true;
}

bool check_nonlocal_dist(const std::vector<Vec3> &pos_trial, const Box &box,
                         const double rh2,
                         const std::optional<double> stair2,
                         const NonlocalBonds &transient_bonds,
                         const NonlocalBonds &permanent_bonds) {
  const unsigned int nbeads = pos_trial.size();


  for (unsigned int i = 0; i < nbeads - 3; i++) {
    for (unsigned int j = i + 3; j < nbeads; j++) {
      double dx = pos_trial[i].x - pos_trial[j].x;
      double dy = pos_trial[i].y - pos_trial[j].y;
      double dz = pos_trial[i].z - pos_trial[j].z;
      box.mindist(dx, dy, dz);

      const double dist2 = dx * dx + dy * dy + dz * dz;

      if (!(dist2 > rh2)) {
        LOG_DEBUG("nonlocal beads " << i << " and " << j
                                    << " overlap with distance "
                                    << std::sqrt(dist2));
        return false;
      }

      const std::tuple<Config, double> p_bond_mask_tuple = permanent_bonds.get_bond_mask(i, j);
      const std::tuple<Config, double> t_bond_mask_tuple = transient_bonds.get_bond_mask(i, j);
      const Config t_bond_mask = std::get<0>(t_bond_mask_tuple);
      const Config p_bond_mask = std::get<0>(p_bond_mask_tuple);

      const double rc2_inner = get_rc2_inner(t_bond_mask_tuple, p_bond_mask_tuple);

#ifdef LOCAL_DEBUG

      if (p_bond_mask)
      {
        std::cout << " Expect permanent bond between " << i << " - " << j << " at rc = " << sqrt(rc2_inner) << " distance = "
            << sqrt(dist2) << std::endl;
      }

#endif

      if (p_bond_mask && !(dist2 < rc2_inner)) {
        LOG_DEBUG("permanent bond between beads "
                  << i << " and " << j << " separated by distance "
                  << std::sqrt(dist2) << " has been broken");
        return false;
      }

      if (stair2 && t_bond_mask && !(dist2 < *stair2)) {
        LOG_DEBUG("transient bond between beads "
                  << i << " and " << j << " separated by distance "
                  << std::sqrt(dist2) << " has exceeded staircase wall");
        return false;
      }
    }
  }

  return true;
}

// calculate the integer of the trial configuration from the post-rotation
// positions of the beads
UpdateConfig config_int(const std::vector<Vec3> &pos_trial, const Box &box,
                        const NonlocalBonds &transient_bonds) {
  const unsigned int nbeads = pos_trial.size();
  UpdateConfig update_config;

  for (unsigned int i = 0; i < nbeads; i++) {
    for (unsigned int j = i + 1; j < nbeads; j++) {


      const std::tuple<Config, double> t_bond_mask_tuple = transient_bonds.get_bond_mask(i, j);
      const Config bond_mask = std::get<0>(t_bond_mask_tuple);
      double rc2 = std::get<1>(t_bond_mask_tuple); // use correct distance

      double dx = pos_trial[i].x - pos_trial[j].x;
      double dy = pos_trial[i].y - pos_trial[j].y;
      double dz = pos_trial[i].z - pos_trial[j].z;
      box.mindist(dx, dy, dz);

      const double dist2 = dx * dx + dy * dy + dz * dz;

      if (dist2 > rc2) {
        continue; 
      }

      // flip bit to form bond
      if (bond_mask != 0) {
        assert(update_config.non_bonded(bond_mask));
        update_config.flip_bond(bond_mask);
        assert(update_config.bonded(bond_mask));
      }
    }
  }

  return update_config;
}

bool accept_move(const std::vector<double> &s_bias, UpdateConfig &orig_config,
                 UpdateConfig &trial_config, Random &mt) {
  double trial_entropy = s_bias[trial_config.config];
  double orig_entropy = s_bias[orig_config.config];
  double dS = trial_entropy - orig_entropy;

  if (dS > 0) {
    std::uniform_real_distribution<> uniform;

    if (uniform(mt) > std::exp(-dS)) {
      return false;
    }
  }

  return true;
}

void crankshaft(std::vector<Vec3> &pos,
                UpdateConfig &update_config,
                const Box &box,
                const double near_min2, const double near_max2,
                const double nnear_min2, const double nnear_max2,
                const double rh2,
                const std::optional<double> stair2,
                const NonlocalBonds &transient_bonds,
                const NonlocalBonds &permanent_bonds, Random &mt,
                const std::vector<double> &s_bias) {
  const unsigned int nbeads = pos.size();

  // trial positions of all beads, since we do not know if the rotated
  // configuration will be accepted
  std::vector<Vec3> pos_trial(pos.size());

  // angle between [0, 360] degrees for rotating the bond between ind and
  // ind - 1
  std::uniform_real_distribution<> theta(0.0, 2 * M_PI);
  // index of a bead lying in the range [1, nbeads - 2]
  std::uniform_int_distribution<> ind(1, nbeads - 2);

  do {
    unsigned int bead_index = ind(mt);
    const double angle = theta(mt);
    rodrigues_rotation(pos, angle, pos_trial, bead_index, box);
  } while (
      !(check_local_dist_if_crankshaft(pos_trial, box, near_min2, near_max2,
                                       nnear_min2, nnear_max2) &&
        check_nonlocal_dist(pos_trial, box, rh2, stair2,
                            transient_bonds, permanent_bonds)));

  UpdateConfig trial_config = config_int(pos_trial, box, transient_bonds);

  if (accept_move(s_bias, update_config, trial_config, mt)) {
    LOG_DEBUG("accept move");
    std::swap(pos, pos_trial);
    std::swap(update_config, trial_config);
  }
}
