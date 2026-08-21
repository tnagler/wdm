// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "methods.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

namespace wdm {

namespace utils {

inline void
remove_incomplete(std::vector<double>& x,
                  std::vector<double>& y,
                  std::vector<double>& w)
{
  size_t complete_count = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    bool row_has_nan = (std::isnan(x[i]) || std::isnan(y[i]));
    if (w.size() > 0)
      row_has_nan = (row_has_nan || std::isnan(w[i]));
    if (!row_has_nan) {
      x[complete_count] = x[i];
      y[complete_count] = y[i];
      if (w.size() > 0)
        w[complete_count] = w[i];
      ++complete_count;
    }
  }

  x.resize(complete_count);
  y.resize(complete_count);
  if (w.size() > 0)
    w.resize(complete_count);
}

inline bool
any_nan(const std::vector<double>& x)
{
  for (size_t i = 0; (i < x.size()); i++) {
    if (std::isnan(x[i]))
      return true;
  }

  return false;
}

inline void
validate_weights(const std::vector<double>& weights)
{
  if (weights.empty())
    return;
  double weight_sum = 0.0;
  for (const auto& weight : weights) {
    if (!std::isfinite(weight) || weight < 0.0)
      throw std::runtime_error("weights must be finite and nonnegative.");
    weight_sum += weight;
  }
  if (!std::isfinite(weight_sum) || weight_sum <= 0.0)
    throw std::runtime_error("weights must have a finite, positive sum.");
}

inline std::string
preproc(std::vector<double>& x,
        std::vector<double>& y,
        std::vector<double>& weights,
        std::string method,
        bool remove_missing)
{
  if (!methods::is_supported(method))
    throw std::runtime_error("method not implemented.");

  if (remove_missing) {
    utils::remove_incomplete(x, y, weights);
    utils::validate_weights(weights);
    if (x.size() < methods::get_min_nobs(method))
      return "return_nan";
  } else {
    std::stringstream msg;
    if (utils::any_nan(x) || utils::any_nan(y) || utils::any_nan(weights)) {
      msg << "there are missing values in the data; "
          << "try remove_missing = TRUE";
    } else {
      utils::validate_weights(weights);
      if (x.size() < methods::get_min_nobs(method))
        msg << "need at least " << methods::get_min_nobs(method)
            << " observations.";
    }
    if (!msg.str().empty())
      throw std::runtime_error(msg.str());
  }

  return "continue";
}

} // end utils

} // end wdm
