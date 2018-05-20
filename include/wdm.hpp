// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "wdm/ktau.hpp"
#include "wdm/hoeffd.hpp"
#include "wdm/prho.hpp"
#include "wdm/srho.hpp"
#include "wdm/bbeta.hpp"
#include "wdm/methods.hpp"
#include "wdm/nan_handling.hpp"

#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/normal.hpp>


//! Weighted dependence measures
namespace wdm {

//! calculates (weighted) dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; see details for possible values.
//! @param weights an optional vector of weights for the data.
//! @param remove_missing if `TRUE`, all observations containing a `nan` are
//!    removed; otherwise throws an error if `nan`s are present.
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
inline double wdm(std::vector<double> x,
                  std::vector<double> y,
                  std::string method,
                  std::vector<double> weights = std::vector<double>(),
                  bool remove_missing = TRUE)
{
    // na handling
    if (utils::preproc(x, y, weights, method, remove_missing) == "return_nan")
        return std::numeric_limits<double>::quiet_NaN();

    if (methods::is_hoeffding(method))
        return hoeffd(x, y, weights);
    if (methods::is_kendall(method))
        return ktau(x, y, weights);
    if (methods::is_pearson(method))
        return prho(x, y, weights);
    if (methods::is_spearman(method))
        return srho(x, y, weights);
    if (methods::is_blomqvist(method))
        return bbeta(x, y, weights);
    throw std::runtime_error("method not implemented.");
}

};


}
