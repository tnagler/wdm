// Copyright © 2025 Thibault Vatter
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "ranks.hpp"
#include "utils.hpp"
#include <tuple>

namespace wdm {
namespace impl {

// Asymptotic standard deviation for xi under null
inline double
xi_std(const std::vector<double>& r,
       const std::vector<double>& l,
       bool y_continuous,
       const std::vector<double>& weights = std::vector<double>())
{
  double n =
    (weights.size() > 0) ? utils::sum(weights) : static_cast<double>(r.size());
  if (y_continuous) {
    return std::sqrt(2.0 / 5.0) / std::sqrt(n);
  } else {
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
}

// Weighted Chatterjee's xi statistic
inline std::tuple<double, double>
cxi(std::vector<double> x,
    std::vector<double> y,
    std::vector<double> weights = std::vector<double>(),
    bool calculate_std = true,
    bool y_continuous = true,
    std::string ties_method = "max")
{
  utils::check_sizes(x, y, weights);

  // Sort x, y, and weights in x order
  utils::sort_all(x, y, weights);

  // Compute ranks of y (ri: number of j such that Y(j) ≤ Y(i)), ties_method
  // "max"
  std::vector<double> r = rank0(y, weights, ties_method);

  // Compute ranks of -y (li: number of j such that Y(j) ≥ Y(i)), ties_method
  // "max"
  std::vector<double> y_neg(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    y_neg[i] = -y[i];
  std::vector<double> l = rank0(y_neg, weights, ties_method);

  // Numerator: sum of weighted absolute differences of consecutive ranks
  double num = 0.0;
  for (size_t i = 1; i < r.size(); ++i)
    num += std::abs(r[i] - r[i - 1]) * (weights.size() > 0 ? weights[i] : 1.0);

  double n =
    (weights.size() > 0) ? utils::sum(weights) : static_cast<double>(x.size());
  double xi;

  if (y_continuous) {
    xi = 1.0 - 3.0 * num / (n * n - 1.0);
  } else {
    // Denominator: sum over (n - li) * li, weighted
    double den = 0.0;
    for (size_t i = 0; i < l.size(); ++i)
      den += (n - l[i]) * l[i] * (weights.size() > 0 ? weights[i] : 1.0);
    xi = 1.0 - n * num / (2.0 * den);
  }

  if (!calculate_std) {
    return std::make_tuple(xi, std::numeric_limits<double>::quiet_NaN());
  } else {
    double std = xi_std(r, l, y_continuous, weights);
    return std::make_tuple(xi, std);
  }
}

} // namespace impl
} // namespace wdm