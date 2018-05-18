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
#include "wdm/indeptest.hpp"
#include "wdm/methods.hpp"

//! Weighted dependence measures
namespace wdm {

//! calculates (weighted) dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"bbeta"`, `"hoeffd"`.
//! @param weights an optional vector of weights for the data.
//! @return the dependence measure
inline double wdm(const std::vector<double>& x,
                  const std::vector<double>& y,
                  std::string method,
                  std::vector<double> weights = std::vector<double>())
{
    using namespace wdm_methods;
    if (is_hoeffding(method))
        return hoeffd(x, y, weights);
    if (is_kendall(method))
        return ktau(x, y, weights);
    if (is_pearson(method))
        return prho(x, y, weights);
    if (is_spearman(method))
        return srho(x, y, weights);
    if (is_blomqvist(method))
        return bbeta(x, y, weights);
    throw std::runtime_error("method not implemented.");
}

}
