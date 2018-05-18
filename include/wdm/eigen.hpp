// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include <Eigen/Dense>
#include "../wdm.hpp"

namespace wdm_eigen {

inline std::vector<double> convert_vec(const Eigen::VectorXd& x)
{
    std::vector<double> xx(x.size());
    Eigen::VectorXd::Map(&xx[0], x.size()) = x;
    return xx;
}

}

namespace wdm {

//! calculates (weighted) dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; see details for possible values. 
//! @param weights an optional vector of weights for the data.
//! 
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation  
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$  
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$  
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$  
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$  
//! 
//! @return the dependence measure
inline double wdm(const Eigen::VectorXd& x,
                  const Eigen::VectorXd& y,
                  std::string method,
                  Eigen::VectorXd weights = Eigen::VectorXd())
{
    using namespace wdm_eigen;
    return wdm(convert_vec(x),
               convert_vec(y),
               method,
               convert_vec(weights));
}

//! calculates a matrix of (weighted) dependence measures.
//! @param x input data.
//! @param method the dependence measure; see details for possible values. 
//! @param weights an optional vector of weights for the data.
//! 
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation  
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$  
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$  
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$  
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$  
//! 
//! @return a matrix of pairwise dependence measures.
inline Eigen::MatrixXd wdm(const Eigen::MatrixXd& x,
                           std::string method,
                           Eigen::VectorXd weights = Eigen::VectorXd())
{
    using namespace wdm_eigen;
    size_t d = x.cols();
    if (d == 1)
        throw std::runtime_error("x must have at least 2 columns.");

    Eigen::MatrixXd ms = Eigen::MatrixXd::Identity(d, d);
    for (size_t i = 0; i < d; i++) {
        for (size_t j = i + 1; j < d; j++) {
            ms(i, j) = wdm(convert_vec(x.col(i)),
                           convert_vec(x.col(j)),
                           method,
                           convert_vec(weights));
            ms(j, i) = ms(i, j);
        }
    }

    return ms;
}

}
