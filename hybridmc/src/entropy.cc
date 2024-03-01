// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// entropy.cc calculates biased entropy, which is approximation to true
// entropy, and uses G-test to detemine if program has converged (which occurs
// when each state has been visited roughly equal number of times); if program
// converges, biased entropy is approximately equal to true entropy

#include "entropy.h"
#include <bit>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>

void init_s(std::vector<double> &s_bias, const unsigned int nbonds) {
  unsigned int nstates = std::pow(2, nbonds);

  for (unsigned int i = 0; i < nstates; i++) {
    s_bias.emplace_back(3 * (nbonds - std::popcount(i)));
  }
}

// calculate the average number of times each configuration is visited (the
// normalization factor)
double compute_norm(const std::vector<uint64_t> &config_count,
                    const unsigned int nstates) {
  double norm = 0.0;
  for (unsigned int i = 0; i < nstates; i++) {
    norm += config_count[i];
  }

  norm /= double(nstates);

  return norm;
}

void compute_entropy(std::vector<uint64_t> &config_count,
                     std::vector<double> &s_bias, const unsigned int nstates,
                     const double pos_scale, const double neg_scale) {
  const double norm = compute_norm(config_count, nstates);

  double s_native = 0.0;
  unsigned int native_ind = nstates - 1;
  double native_count = config_count[native_ind];

  // if the native state is visited,
  if (native_count > 0) {
    // calculate the entropy of the native state
    s_native = s_bias[native_ind] + pos_scale * std::log(native_count / norm);
  } else {
    s_native = s_bias[native_ind] - neg_scale * std::log(norm);
  }

  for (unsigned int i = 0; i < nstates; i++) {
    double s_i = s_bias[i];

    // if a configuration is visited,
    if (config_count[i] > 0) {
      // calculate the entropy of the configuration
      s_i = s_bias[i] + pos_scale * std::log(config_count[i] / norm);
    } else {
      // lower the entropy to make unobserved states more likely
      s_i = s_bias[i] - neg_scale * std::log(norm);
    }

    s_bias[i] = s_i - s_native;
    std::cout << "entropy of " << i << " is " << s_bias[i] << " with count "
              << config_count[i] << std::endl;
  }
}

// how to calculate G value: https://en.wikipedia.org/wiki/G-test
double compute_g_val(std::vector<uint64_t> &config_count,
                     const unsigned int nstates) {
  // test statistic
  double norm = compute_norm(config_count, nstates);

  double g_val = 0.0;

  // observed value of test statistic for set of entropies
  for (unsigned int i = 0; i < nstates; i++) {
    if (config_count[i] > 0) {
      g_val += config_count[i] * std::log(config_count[i] / norm);
    }
  }
  g_val *= 2.0;

  return g_val;
}

// use G-test to test for uniformity
bool g_test(std::vector<uint64_t> &config_count, const unsigned int nstates,
            const double sig_level) {
  unsigned int dof = nstates - 1;
  double g_val = compute_g_val(config_count, nstates);

  // find the probability of observing a statistic (g_val) as extreme as the
  // test statistic (norm) using an upper-tailed test
  boost::math::chi_squared chi_dist(dof);
  double p_val = 1.0 - boost::math::cdf(chi_dist, g_val);

  // critical value of G, the threshold beyond which lies the region of a
  // distribution where results are statistically significant
  double g_crit =
      boost::math::quantile(boost::math::complement(chi_dist, sig_level));

  if (!(g_val < g_crit)) {
    std::cout << "Rejected uniform test: g_val " << g_val << ", g_crit " << g_crit
              << ", p_val " << p_val << std::endl;
    return false;
  }

  std::cout << "Acceptable uniform test: g_val " << g_val << ", g_crit " << g_crit
            << ", p_val " << p_val << std::endl;

  return true;
}
