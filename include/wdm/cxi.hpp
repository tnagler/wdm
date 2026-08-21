// Copyright © 2025 Thibault Vatter
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "ranks.hpp"
#include "utils.hpp"
#include <memory>
#include <tuple>

namespace wdm {
namespace impl {

//! Seeds used to break predictor ties when the caller supplies none. A
//! dependence measure has to be a function of its arguments, so the default is
//! a constant rather than a draw from `std::random_device`.
inline std::vector<int>
default_tie_seeds()
{
  return { 1, 2, 3, 4, 5 };
}

//! Sort observations by the predictor and break predictor ties uniformly at
//! random, independently of the response.
inline void
sort_chatterjee_observations(std::vector<double>& x,
                             std::vector<double>& y,
                             std::vector<double>& weights,
                             const std::vector<int>& seeds)
{
  std::vector<size_t> order = utils::get_order(x);
  std::unique_ptr<random::RandomGenerator> tie_generator;
  for (size_t begin = 0, end; begin < order.size(); begin = end) {
    end = begin + 1;
    while (end < order.size() && x[order[end]] == x[order[begin]])
      ++end;
    if (end - begin > 1) {
      if (!tie_generator)
        tie_generator.reset(new random::RandomGenerator(seeds));
      std::vector<size_t> tied_order(order.begin() + begin,
                                     order.begin() + end);
      random::shuffle(tied_order, *tie_generator);
      std::copy(tied_order.begin(), tied_order.end(), order.begin() + begin);
    }
  }

  std::vector<double> sorted_x(x.size()), sorted_y(y.size()),
    sorted_weights(weights.size());
  for (size_t i = 0; i < order.size(); ++i) {
    sorted_x[i] = x[order[i]];
    sorted_y[i] = y[order[i]];
    sorted_weights[i] = weights[order[i]];
  }
  x = sorted_x;
  y = sorted_y;
  weights = sorted_weights;
}

// Conditional null mean and standard deviation for a continuous response.
inline std::tuple<double, double>
xi_continuous_inference(const std::vector<double>& probabilities)
{
  double edge_weight_sum = 0.0;
  double null_numerator_mean = 0.0;
  double squared_edge_weight_sum = 0.0;
  double adjacent_edge_product_sum = 0.0;
  double squared_probability_sum = 0.0;
  double edge_node_product_sum = 0.0;

  for (size_t i = 0; i < probabilities.size(); ++i)
    squared_probability_sum += probabilities[i] * probabilities[i];

  for (size_t i = 0; i + 1 < probabilities.size(); ++i) {
    edge_weight_sum += probabilities[i];
    null_numerator_mean +=
      probabilities[i] *
      (1.0 / 3.0 + (probabilities[i] + probabilities[i + 1]) / 6.0);
    squared_edge_weight_sum += probabilities[i] * probabilities[i];
    edge_node_product_sum +=
      probabilities[i] * (probabilities[i] + probabilities[i + 1]);
    if (i + 2 < probabilities.size())
      adjacent_edge_product_sum += probabilities[i] * probabilities[i + 1];
  }

  double null_numerator_variance = squared_edge_weight_sum / 18.0;
  null_numerator_variance += adjacent_edge_product_sum / 90.0;
  null_numerator_variance +=
    edge_weight_sum * edge_weight_sum * squared_probability_sum / 45.0;
  null_numerator_variance -= edge_weight_sum * edge_node_product_sum / 45.0;
  if (!std::isfinite(null_numerator_variance) || null_numerator_variance <= 0.0)
    throw std::runtime_error(
      "cannot compute the null variance of Chatterjee's xi.");

  return std::make_tuple(3.0 * std::sqrt(null_numerator_variance),
                         1.0 - 3.0 * null_numerator_mean);
}

// Asymptotic standard deviation for xi with a tied response.
inline double
xi_std(const std::vector<double>& r,
       const std::vector<double>& l,
       const std::vector<double>& weights = std::vector<double>())
{
  double n =
    (weights.size() > 0) ? utils::sum(weights) : static_cast<double>(r.size());

  // Weighted version
  std::vector<double> i(r.size());
  for (size_t k = 0; k < r.size(); ++k)
    i[k] = k + 1;

  // Sort r and weights together
  std::vector<size_t> order = utils::get_order(r);
  std::vector<double> u(r.size()), w(r.size());
  for (size_t k = 0; k < r.size(); ++k) {
    u[k] = r[order[k]];
    w[k] = (weights.size() > 0) ? weights[order[k]] : 1.0;
  }

  // Weighted cumulative sum
  std::vector<double> v(r.size());
  v[0] = u[0] * w[0];
  for (size_t k = 1; k < r.size(); ++k)
    v[k] = v[k - 1] + u[k] * w[k];

  double an = 0, bn = 0, cn = 0, dn = 0;
  for (size_t k = 0; k < r.size(); ++k) {
    an += (2 * n - 2 * i[k] + 1) * u[k] * u[k] * w[k];
    cn += (2 * n - 2 * i[k] + 1) * u[k] * w[k];
    dn += l[k] * (n - l[k]) * ((weights.size() > 0) ? weights[k] : 1.0);
  }
  an /= std::pow(n, 4);
  cn /= std::pow(n, 3);
  dn /= std::pow(n, 3);

  for (size_t k = 0; k < r.size(); ++k) {
    double temp = v[k] + (n - i[k]) * u[k] * w[k];
    bn += temp * temp;
  }
  bn /= std::pow(n, 5);

  double tau2 = (an - 2 * bn + cn * cn) / (dn * dn);
  return std::sqrt(tau2) / std::sqrt(n);
}

//! Weighted Chatterjee's xi statistic and conditional null inference.
//! @param x predictor values.
//! @param y response values.
//! @param weights optional case weights, normalized internally.
//! @param calculate_std whether to calculate analytic null inference.
//! @param ties_method rank convention for tied responses.
//! @param seeds optional seeds for random predictor-tie breaking.
//! @param y_continuous whether the response distribution is known to be
//!   continuous. Observed response ties always select tied-response inference.
//! @return `(estimate, standard_error, null_mean, inference_estimate)`. The
//!   inferential values are `NaN` when `calculate_std` is false.
//! @details Weights must be finite, nonnegative, and have a positive sum.
//!   Analytic inference with unequal weights assumes a continuous response and
//!   weights that are fixed or depend only on `x`. With unequal weights, it is
//!   unavailable for a discrete or tied response. For a continuous response,
//!   `inference_estimate` is \f$1 - 3 A\f$; otherwise it is the reported
//!   denominator-corrected estimate.
inline std::tuple<double, double, double, double>
cxi(std::vector<double> x,
    std::vector<double> y,
    std::vector<double> weights = std::vector<double>(),
    bool calculate_std = true,
    std::string ties_method = "max",
    std::vector<int> seeds = std::vector<int>(),
    bool y_continuous = true)
{
  utils::check_sizes(x, y, weights);

  if (weights.size() == 0)
    weights = std::vector<double>(x.size(), 1.0);

  utils::validate_weights(weights);
  double weight_sum = utils::sum(weights);

  // Zero-mass observations are absent from the weighted empirical measure and
  // must not create additional edges in predictor order.
  for (size_t i = weights.size(); i-- > 0;) {
    if (weights[i] == 0.0) {
      x.erase(x.begin() + i);
      y.erase(y.begin() + i);
      weights.erase(weights.begin() + i);
    }
  }

  // Sort in x order and break x ties uniformly without consulting y. An empty
  // seed vector would draw from std::random_device, making the estimate differ
  // between calls on the same data; pass seeds to vary the tie ordering.
  sort_chatterjee_observations(
    x, y, weights, seeds.empty() ? default_tie_seeds() : seeds);

  std::vector<double> probabilities = weights;
  for (auto& probability : probabilities)
    probability /= weight_sum;
  bool weights_are_unequal = false;
  for (size_t i = 1; i < probabilities.size(); ++i)
    weights_are_unequal =
      weights_are_unequal || probabilities[i] != probabilities[0];

  std::vector<double> ordered_response = y;
  std::sort(ordered_response.begin(), ordered_response.end());
  if (ordered_response.front() == ordered_response.back())
    throw std::runtime_error(
      "Chatterjee's xi is undefined for a constant response.");
  bool response_has_ties =
    std::adjacent_find(ordered_response.begin(), ordered_response.end()) !=
    ordered_response.end();

  // Weighted empirical distribution at each response.
  std::vector<double> r = rank0(y, probabilities, ties_method);

  // Weighted empirical survival function at each response.
  std::vector<double> y_neg(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    y_neg[i] = -y[i];
  std::vector<double> l = rank0(y_neg, probabilities, ties_method);

  // Numerator: base-point weight on edge (i, i + 1).
  double num = 0.0;
  for (size_t i = 0; i + 1 < r.size(); ++i)
    num += probabilities[i] * std::abs(r[i + 1] - r[i]);

  // General weighted-rank denominator, valid for continuous and tied responses.
  double den = 0.0;
  for (size_t i = 0; i < l.size(); ++i)
    den += 2.0 * probabilities[i] * l[i] * (1.0 - l[i]);
  if (!std::isfinite(den) || den <= 0.0)
    throw std::runtime_error(
      "Chatterjee's xi is undefined for a constant response.");

  double xi = 1.0 - num / den;

  if (!calculate_std) {
    return std::make_tuple(xi,
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN());
  } else if (y_continuous && !response_has_ties) {
    auto inference = xi_continuous_inference(probabilities);
    return std::make_tuple(
      xi, std::get<0>(inference), std::get<1>(inference), 1.0 - 3.0 * num);
  } else {
    if (weights_are_unequal)
      throw std::runtime_error(
        "analytic Chatterjee inference is unavailable for an unequally "
        "weighted, discrete or tied response.");
    std::vector<double> raw_r = rank0(y, {}, ties_method);
    std::vector<double> raw_l = rank0(y_neg, {}, ties_method);
    return std::make_tuple(xi, xi_std(raw_r, raw_l), 0.0, xi);
  }
}

} // namespace impl
} // namespace wdm
