// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "utils.hpp"
#include "ktau.hpp"
#include "hoeffd.hpp"
#include "prho.hpp"
#include "srho.hpp"
#include "bbeta.hpp"
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/normal.hpp>
#include <random>
#include <limits>
#include "methods.hpp"

namespace wdm_interp {
    
inline double linear_interp(const double& x,
                     const std::vector<double>& grid,
                     const std::vector<double>& values)
{
    // find upper end point of interval
    size_t i = 1;
    while(x > grid[i]) i++;
    
    // linear interpolation
    double w = (x - grid[i - 1]) / (grid[i] - grid[i - 1]);
    return w * values[i - 1] + (1 - w) * values[i];
}

}


namespace wdm {

//! calculates the (approximate) asymptotic distribution function of Hoeffding's
//! B (as in Blum, Kiefer, and Rosenblatt) under the null hypothesis of
//! independence.
//! @param B sample estimate of Hoeffding's B.
//! @param n the sample size.
inline double phoeffb(double B, size_t n) {
    B *= 0.5 * std::pow(boost::math::constants::pi<double>(), 4) * (n - 1);

    // obtain approximate p values by interpolation of tabulated values
    using namespace boost::math;
    double p;
    if ((B <= 1.1) | (B >= 8.5)) {
        p = std::min(1.0, std::exp(0.3885037 - 1.164879 * B));
        p = std::max(1e-12, p);
    } else {
        std::vector<double> grid{
            1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6,
            1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2,
            2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75,
            2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35,
            3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95,
            4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55,
            4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5, 5.5, 6, 6.5, 7,
            7.5, 8, 8.5
        };
        std::vector<double> vals{
            0.5297, 0.4918, 0.4565, 0.4236, 0.3930, 0.3648, 0.3387, 0.3146,
            0.2924, 0.2719, 0.2530, 0.2355, 0.2194, 0.2045, 0.1908, 0.1781,
            0.1663, 0.1554, 0.1453, 0.1359, 0.1273, 0.1192, 0.1117, 0.1047,
            0.0982, 0.0921, 0.0864, 0.0812, 0.0762, 0.0716, 0.0673, 0.0633,
            0.0595, 0.0560, 0.0527, 0.0496, 0.0467, 0.0440, 0.0414, 0.0390,
            0.0368, 0.0347, 0.0327, 0.0308, 0.0291, 0.0274, 0.0259, 0.0244,
            0.0230, 0.0217, 0.0205, 0.0194, 0.0183, 0.0173, 0.0163, 0.0154,
            0.0145, 0.0137, 0.0130, 0.0123, 0.0116, 0.0110, 0.0104, 0.0098,
            0.0093, 0.0087, 0.0083, 0.0078, 0.0074, 0.0070, 0.0066, 0.0063,
            0.0059, 0.0056, 0.0053, 0.0050, 0.0047, 0.0045, 0.0042, 0.00025,
            0.00014, 0.0008, 0.0005, 0.0003, 0.0002, 0.0001
        };
        p = wdm_interp::linear_interp(B, grid, vals);
    }

    return p;
}

//! calculates the test statistic for indpendence tests
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
inline double calculate_test_stat(
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::string method,
        std::vector<double> weights = std::vector<double>())
{
    // determine effective sample size
    double n_eff = wdm_utils::effective_sample_size(x.size(), weights);

    // calculate test statistic
    using namespace wdm_methods;
    double stat;
    if (is_hoeffding(method)) {
        stat = hoeffd(x, y, weights) / 30.0 + 1.0 / (36.0 * n_eff);
    } else if (is_kendall(method)) {
        stat = ktau(x, y, weights);
        stat *= std::sqrt(9 * n_eff / 4);
    } else if (is_pearson(method)) {
        stat = boost::math::atanh(prho(x, y, weights));
        stat *= std::sqrt(n_eff - 3);
    } else if (is_spearman(method)) {
        stat = boost::math::atanh(srho(x, y, weights));
        stat *= std::sqrt((n_eff - 3) / 1.06);
    }  else if (is_blomqvist(method)) {
        stat = bbeta(x, y, weights);
        stat *= std::sqrt(n_eff);
    } else {
        throw std::runtime_error("method not implemented.");
    }

    return stat;
}

//! calculates the asymptotic p-value.
//! @param stat value of the test statistic.
//! @param method the dependence measure; see details for possible values. 
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation  
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$  
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$  
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$  
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$  
//! 
//! @param n_eff effective sample size; only used for method `"hoeffd"`.
inline double calculate_asymptotic_p_val(double stat,
                                         std::string method,
                                         double n_eff = 0.0)
{
    double p_val;
    if (method == "hoeffd") {
        if (n_eff == 0.0)
            throw std::runtime_error("must provide n_eff for method 'hoeffd'.");
        p_val = phoeffb(stat, n_eff);
    } else {
        boost::math::normal norm_dist(0, 1);
        p_val = 2 * boost::math::cdf(norm_dist, -std::abs(stat));
    }

    return p_val;
}

//! calculates asymptotic p-values of independence tests based on (weighted)
//! dependence measures.
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
//! @return the p-value of the independence test.
inline double indeptest(
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::string method,
        std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);
    double stat = calculate_test_stat(x, y, method, weights);
    double n_eff = wdm_utils::effective_sample_size(x.size(), weights);

    return calculate_asymptotic_p_val(stat, method, n_eff);
}

}
