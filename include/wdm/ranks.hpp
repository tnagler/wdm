// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "nan_handling.hpp"
#include "random.hpp"
#include "utils.hpp"

#include <algorithm>
#include <numeric>

namespace wdm {

namespace impl {

//! draws a key for every observation, the order in which `rank(x, weights,
//! "random", seeds)` breaks ties.
//!
//! Within a group of tied values, the observations receive the group's ranks
//! in increasing order of their keys. The keys are a uniformly random
//! permutation, so every group is ordered uniformly at random, independently
//! of the others; and an observation's key depends only on `n` and `seeds`, so
//! a group's order depends only on which observations are in it. Values that
//! gain or lose a tie leave every other group's order as it was.
//! @param n number of observations.
//! @param seeds seeds of the random number generator, as for `rank()`.
//! @return a permutation of `0, ..., n - 1`: observation `i`'s key.
inline std::vector<size_t>
tie_keys(size_t n, std::vector<int> seeds = std::vector<int>())
{
  std::vector<size_t> keys(n);
  std::iota(keys.begin(), keys.end(), 0);
  if (n > 1) {
    random::RandomGenerator generator(seeds);
    random::shuffle(keys, generator);
  }
  return keys;
}

//! computes ranks.
//! @param x input vector.
//! @param weights (optional), weights for each observation.
//! @param ties_method `"min"` (default) assigns all tied values the minimum
//!   score; `"average"` assigns the average score, `"first"` ranks them in
//!   order of occurance, `"random"` in a uniformly random order, drawn as the
//!   order of the observations' keys (see `tie_keys()`).
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return a vector containing the ranks of each element in `x`.
inline std::vector<double>
rank(std::vector<double> x,
     std::vector<double> weights = std::vector<double>(),
     std::string ties_method = "min",
     std::vector<int> seeds = std::vector<int>())
{
  if ((ties_method != "min") && (ties_method != "average") &&
      (ties_method != "first") && (ties_method != "random"))
    throw std::runtime_error(
      "ties method must be one of 'min', 'average', 'first', 'random'.");

  // set default weights if necessary
  size_t n = x.size();
  if (weights.size() == 0)
    weights = std::vector<double>(n, 1.0);

  if (weights.size() != n) {
    throw std::runtime_error("weights and data must have same size.");
  }

  // NaN-handling
  std::vector<double> nans;
  if (utils::any_nan(x)) {
    nans.resize(n, 0);
    for (size_t i = 0; i < n; i++) {
      if (std::isnan(x[i])) {
        x[i] = std::numeric_limits<double>::max();
        nans[i] = 1;
        weights[i] = 0;
      }
    }
  }

  double w_mean =
    utils::sum(weights) / static_cast<double>(n - utils::sum(nans));
  for (auto& w : weights) {
    w = w / w_mean;
  }

  // permutation that brings 'x' in ascending order; under "random", ties are
  // ordered by the observations' keys, which makes the order total, so that
  // no sort algorithm can arrange it differently
  std::vector<size_t> perm;
  if (ties_method == "random") {
    const std::vector<size_t> keys = tie_keys(n, seeds);
    perm.resize(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b) {
      return (x[a] < x[b]) || ((x[a] == x[b]) && (keys[a] < keys[b]));
    });
  } else {
    perm = utils::get_order(x);
  }

  double w_acc = 0.0, w_batch;
  for (size_t i = 0, reps; i < n; i += reps) {
    // find replications
    reps = 0;
    w_batch = 0.0;
    while ((i + reps < n) && (x[perm[i]] == x[perm[i + reps]]))
      w_batch += weights[perm[i + reps++]];

    // assign min rank
    for (size_t k = 0; k < reps; ++k)
      x[perm[i + k]] = w_acc + weights[perm[i]];

    if (reps > 1) {
      if ((ties_method == "first") || (ties_method == "random")) {
        // break ties by assigning the cumulative weights, in order of
        // appearance ("first") or of the keys ("random"), which is the order
        // of `perm` either way
        double ww = 0.0;
        for (size_t k = 0; k < reps; ++k) {
          ww += weights[perm[i + k]];
          x[perm[i + k]] = w_acc + ww;
        }
      } else if (ties_method == "average") {
        // assign average rank to tied values
        for (size_t k = 0; k < reps; ++k)
          x[perm[i + k]] += (w_batch - weights[perm[i]]) / 2;
      }
    }

    // accumulate weights for current batch
    w_acc += w_batch;
  }

  if (nans.size() == n) {
    for (size_t i = 0; i < x.size(); i++) {
      if (nans[i]) {
        x[i] = NAN;
      }
    }
  }

  return x;
}

//! computes ranks (such that smallest element has rank 0), assigning average
//! ranks for ties.
//! @param x input vector.
//! @param ties_method `"min"` (default) assigns all tied values the minimum
//!   score; `"average"` assigns the average score; `"max"` assigns the
//!   maximum score, so that a rank is the total weight of the observations
//!   that are less than or equal to the corresponding value.
//! @param weights (optional), weights for each observation.
//! @return a vector containing the ranks of each element in `x`.
inline std::vector<double>
rank0(std::vector<double> x,
      std::vector<double> weights = std::vector<double>(),
      std::string ties_method = "min")
{
  if ((ties_method != "min") && (ties_method != "average") &&
      (ties_method != "max"))
    throw std::runtime_error(
      "ties_method must be either 'min', 'average', or 'max'.");

  // set default weights if necessary
  size_t n = x.size();
  if (weights.size() == 0)
    weights = std::vector<double>(n, 1.0);

  // permutation that brings 'x' in ascending order
  std::vector<size_t> perm = utils::get_order(x);

  double w_acc = 0.0, w_batch;
  for (size_t i = 0, reps; i < n; i += reps) {
    // find replications
    reps = 0;
    w_batch = 0.0;
    while ((i + reps < n) && (x[perm[i]] == x[perm[i + reps]]))
      w_batch += weights[perm[i + reps++]];

    // assign min rank
    for (size_t k = 0; k < reps; ++k)
      x[perm[i + k]] = w_acc;

    // accumulate weights for current batch
    w_acc += w_batch;

    // assign average rank to tied values
    if ((ties_method == "average") && (reps > 1)) {
      std::vector<double> ww(reps);
      for (size_t k = 0; k < reps; ++k)
        ww[k] = weights[perm[i + k]];
      double offset = utils::perm_sum(ww, 2) / w_batch;
      for (size_t k = 0; k < reps; ++k)
        x[perm[i + k]] += offset;
    } else if (ties_method == "max") {
      // w_acc now holds the weight of everything up to and including the batch
      for (size_t k = 0; k < reps; ++k)
        x[perm[i + k]] = w_acc;
    }
  }

  return x;
}

//! computes the bivariate rank of a pair of vectors: for each observation,
//! the (weighted) number of observations whose `x` and `y` are both strictly
//! smaller.
//! @param x first input vector.
//! @param y second input vecotr.
//! @param weights (optional), weights for each observation.
inline std::vector<double>
bivariate_rank(const std::vector<double>& x,
               const std::vector<double>& y,
               std::vector<double> weights = std::vector<double>())
{
  utils::check_sizes(x, y, weights);
  size_t n = x.size();
  if (weights.size() == 0)
    weights = std::vector<double>(n, 1.0);

  // Visit the observations by increasing x, and those with equal x by
  // decreasing y: an observation visited earlier then has a strictly smaller
  // x whenever it has a strictly smaller y.
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
    return (x[i] < x[j]) || ((x[i] == x[j]) && (y[i] > y[j]));
  });

  // The distinct values of y in increasing order, and the weight seen so far
  // at each.
  std::vector<double> levels = y;
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  utils::FenwickTree seen(levels.size());

  std::vector<double> counts(n);
  for (size_t i : order) {
    // Which distinct y-value does this observation have?
    size_t level = static_cast<size_t>(
      std::lower_bound(levels.begin(), levels.end(), y[i]) - levels.begin());
    // How much weight has already been seen at smaller y-values?
    counts[i] = seen.prefix_sum(level);
    // Record this observation at its y-value.
    seen.add(level, weights[i]);
  }

  return counts;
}

//! computes the (weighted) median of a vector.
//! @param x the input vector.
inline double
median(const std::vector<double>& x,
       std::vector<double> weights = std::vector<double>())
{
  utils::check_sizes(x, x, weights);
  size_t n = x.size();

  // sort x and weights in x order
  auto perm = utils::get_order(x);
  auto xx = x;
  auto w = weights;
  for (size_t i = 0; i < n; i++) {
    xx[i] = x[perm[i]];
    if (w.size() > 0)
      w[i] = weights[perm[i]];
  }

  // compute weighted ranks and the "average rank" (corresponds to the
  // median)
  auto ranks = rank0(xx, w, "average");
  if (weights.size() == 0)
    weights = std::vector<double>(n, 1.0);
  double rank_avrg = utils::perm_sum(weights, 2) / utils::sum(weights);

  // weighted median splits data below and above rank_avrg
  size_t i = 0;
  while (ranks[i] < rank_avrg)
    i++;
  if (ranks[i] == rank_avrg)
    return xx[i];
  else
    return 0.5 * (xx[i - 1] + xx[i]);
}
}
}
