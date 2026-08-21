// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "../wdm.hpp"
#include <Eigen/Dense>

namespace wdm {

namespace utils {

inline std::vector<double>
convert_vec(const Eigen::VectorXd& x)
{
  std::vector<double> xx(x.size());
  if (x.size() > 0)
    Eigen::VectorXd::Map(&xx[0], x.size()) = x;
  return xx;
}
}

//! Calculates a dependence measure for two Eigen vectors.
//! @param x, y input vectors of equal length.
//! @param method any method name accepted by the `std::vector` overload.
//! @param weights optional case weights subject to the same validation as the
//!   `std::vector` overload.
//! @param remove_missing whether to remove rows containing a `NaN`.
//! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
//! @return The result of the corresponding `std::vector` overload.
//! @throws std::runtime_error under the same conditions as that overload.
inline double
wdm(const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    std::string method,
    Eigen::VectorXd weights = Eigen::VectorXd(),
    bool remove_missing = true,
    std::vector<int> seeds = std::vector<int>())
{
  return wdm(utils::convert_vec(x),
             utils::convert_vec(y),
             method,
             utils::convert_vec(weights),
             remove_missing,
             seeds);
}

//! Calculates all pairwise dependence measures between matrix columns.
//! @param x matrix with observations in rows and variables in columns.
//! @param method any method name accepted by the vector overload.
//! @param weights optional case weights for the rows.
//! @param remove_missing whether each pair should remove rows containing a
//!   `NaN`.
//! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
//! @return A matrix of pairwise dependence measures. Chatterjee matrices are
//!   generally asymmetric; all other supported measures produce symmetric
//!   matrices.
//! @throws std::runtime_error if `x` has exactly one column or a pairwise call
//!   rejects its inputs.
inline Eigen::MatrixXd
wdm(const Eigen::MatrixXd& x,
    std::string method,
    Eigen::VectorXd weights = Eigen::VectorXd(),
    bool remove_missing = true,
    std::vector<int> seeds = std::vector<int>())
{
  size_t d = x.cols();
  if (d == 1)
    throw std::runtime_error("x must have at least 2 columns.");

  Eigen::MatrixXd ms = Eigen::MatrixXd::Identity(d, d);
  for (size_t i = 0; i < d; i++) {
    for (size_t j = i + 1; j < d; j++) {
      ms(i, j) = wdm(utils::convert_vec(x.col(i)),
                     utils::convert_vec(x.col(j)),
                     method,
                     utils::convert_vec(weights),
                     remove_missing,
                     seeds);
      if (methods::is_chatterjee(method)) {
        ms(j, i) = wdm(utils::convert_vec(x.col(j)),
                       utils::convert_vec(x.col(i)),
                       method,
                       utils::convert_vec(weights),
                       remove_missing,
                       seeds);
      } else {
        ms(j, i) = ms(i, j);
      }
    }
  }

  return ms;
}

}
