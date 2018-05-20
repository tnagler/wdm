// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include <limits>

namespace wdm {

namespace utils {

inline void remove_incomplete(std::vector<double>& x,
                              std::vector<double>& y,
                              std::vector<double>& w)
{
    // if observation conatins nan, set all to nan
    for (size_t i = 0; i < x.size(); i++) {
        bool row_has_nan = (std::isnan(x[i]) | std::isnan(y[i]));
        if (w.size() > 0)
            row_has_nan = (row_has_nan |  std::isnan(w[i]));
        if (row_has_nan) {
            x[i] = std::numeric_limits<double>::quiet_NaN();
            y[i] = std::numeric_limits<double>::quiet_NaN();
            if (w.size() > 0)
                w[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }

    // remove all nan observations
    auto is_nan_pred = [] (const double& val) {
        return std::isnan(val);
    };
    x.erase(std::remove_if(x.begin(), x.end(), is_nan_pred), x.end());
    y.erase(std::remove_if(y.begin(), y.end(), is_nan_pred), y.end());
    w.erase(std::remove_if(w.begin(), w.end(), is_nan_pred), w.end());
}

inline bool any_nan(const std::vector<double>& x) {
    for (size_t i = 0; (i < x.size()); i++) {
        if (std::isnan(x[i]))
            return true;
    }

    return false;
}

inline std::string preproc(std::vector<double>& x,
                           std::vector<double>& y,
                           std::vector<double>& weights,
                           std::string method,
                           bool remove_missing)
{
    size_t min_nobs = (method == "hoeffding") ? 5 : 2;
    if (remove_missing) {
        utils::remove_incomplete(x, y, weights);
        if (x.size() < min_nobs)
            return "return_nan";
    } else {
        std::stringstream msg;
        if (utils::any_nan(x) | utils::any_nan(y) | utils::any_nan(weights)) {
            msg << "there are missing values in the data; " <<
                   "try remove_missing = TRUE";
        } else if (x.size() < min_nobs) {
            msg << "need at least " << min_nobs << "observations.";
        }
        if (!msg.str().empty())
            throw std::runtime_error(msg.str());
    }

    return "continue";
}

} // end utils

} // end wdm
