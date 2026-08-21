// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

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

//! calculates (weighted) dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; see details for possible values.
//! @param weights an optional vector of weights for the data.
//! @param remove_missing if `true`, all observations containing a `nan` are
//!    removed; otherwise throws an error if `nan`s are present.
//! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$
//!   - `"chatterjee"`, `"cxi"`, `"xi"`: Chatterjee's \f$ \xi \f$
//!
//! @return the dependence measure
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

//! calculates a matrix of (weighted) dependence measures.
//! @param x input data.
//! @param method the dependence measure; see details for possible values.
//! @param weights an optional vector of weights for the data.
//! @param remove_missing if `true`, all observations containing a `nan` are
//!    removed; otherwise throws an error if `nan`s are present.
//! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$
//!   - `"chatterjee"`, `"cxi"`, `"xi"`: Chatterjee's \f$ \xi \f$
//!
//! @return a matrix of pairwise dependence measures. Chatterjee matrices are
//!   generally asymmetric; all other supported measures produce symmetric
//!   matrices.
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
