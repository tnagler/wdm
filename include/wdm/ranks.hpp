// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "nan_handling.hpp"
#include "utils.hpp"
#include <random>

namespace wdm {

namespace impl {

//! computes ranks (such that smallest element has rank 0).
//! @param x input vector.
//! @param ties_method `"min"` (default) assigns all tied values the minimum
//!   score; `"average"` assigns the average score, `"first"` ranks them in
//!   order of occurance, `"random"` randomizes.
//! @param weights (optional), weights for each observation.
//! @return a vector containing the ranks of each element in `x`.
inline std::vector<double>
rank(std::vector<double> x,
     std::vector<double> weights = std::vector<double>(),
     std::string ties_method = "min")
{
    if ((ties_method != "min") && (ties_method != "average") &&
        (ties_method != "first") && (ties_method != "random"))
        throw std::runtime_error(
          "ties method must be one of 'min', 'average', 'first', 'random'.");

    // set default weights if necessary
    size_t n = x.size();
    if (weights.size() == 0)
        weights = std::vector<double>(n, 1.0);

    if (weights.size() != n) {
        throw std::runtime_error("weights and data must have same size.");
    }
    
    // NaN-handling
    std::vector<bool> nans;
    if (utils::any_nan(x)) {
        nans.resize(n, false);
        for (size_t i = 0; i < n; i++) {
            if (std::isnan(x[i])) {
                x[i] = std::numeric_limits<double>::max();
                nans[i] = true;
                weights[i] = 0;
            }
        }
    }

    // permutation that brings 'x' in ascending order
    std::vector<size_t> perm = utils::get_order(x);

    double w_acc = 0.0, w_batch;
    for (size_t i = 0, reps; i < n; i += reps) {
        // find replications
        reps = 0;
        w_batch = 0.0;
        while ((i + reps < n) && (x[perm[i]] == x[perm[i + reps]]))
            w_batch += weights[perm[i + reps++]];

        // assign min rank
        for (size_t k = 0; k < reps; ++k)
            x[perm[i + k]] = w_acc;

        // accumulate weights for current batch
        w_acc += w_batch;

        if (reps <= 1)
            continue;

        if (ties_method == "first") {
            // assign weighted ranks in order of appearance
            double ww = 0;
            for (size_t k = 0; k < reps; ++k) {
                ww += weights[perm[i + k]];
                x[perm[i + k]] += ww;
            }
        } else if (ties_method == "random") {
            // assign weighted ranks in random order
            std::random_device rd;
            std::default_random_engine gen(rd());
            std::vector<size_t> rvals(reps);
            std::iota(rvals.begin(), rvals.end(), 0); // 0, 1, 2, ...
            std::shuffle(rvals.begin(), rvals.end(), gen);

            double ww = 0;
            for (size_t k = 0; k < reps; ++k) {
                x[perm[i + rvals[k]]] += ww;
                ww += weights[perm[i + rvals[k]]];
            }
        } else if (ties_method == "average") {
            // assign average rank to tied values
            std::vector<double> ww(reps);
            for (size_t k = 0; k < reps; ++k)
                ww[k] = weights[perm[i + k]];
            for (size_t k = 0; k < reps; ++k)
                x[perm[i + k]] += utils::perm_sum(ww, 2) / w_batch;
        }
    }

    if (nans.size() == n) {
        for (size_t i = 0; i < x.size(); i++) {
            if (nans[i]) {
                x[i] = NAN;
            }
        }
    }

    return x;
}

//! computes the bivariate rank of a pair of vectors (starting at 0).
//! @param x first input vector.
//! @param y second input vecotr.
//! @param weights (optional), weights for each observation.
inline std::vector<double>
bivariate_rank(std::vector<double> x,
               std::vector<double> y,
               std::vector<double> weights = std::vector<double>())
{
    utils::check_sizes(x, y, weights);

    // get inverse of permutation that brings x in ascending order
    std::vector<size_t> perm_x = utils::get_order(x);
    perm_x = utils::invert_permutation(perm_x);

    // sort x, y, and weights according to x, breaking ties with y
    utils::sort_all(x, y, weights);

    // get inverse of permutation that brings y in descending order
    std::vector<size_t> perm_y = utils::get_order(y, false);
    perm_y = utils::invert_permutation(perm_y);

    // sort y in descending order counting inversions
    std::vector<double> counts(y.size(), 0.0);
    utils::merge_sort_count_per_element(y, weights, counts);

    // bring counts back in original order
    std::vector<double> counts_tmp = counts;
    for (size_t i = 0; i < counts.size(); i++)
        counts[i] = counts_tmp[perm_y[perm_x[i]]];

    return counts;
}

//! computes the (weighted) median of a vector.
//! @param x the input vector.
inline double
median(const std::vector<double>& x,
       std::vector<double> weights = std::vector<double>())
{
    utils::check_sizes(x, x, weights);
    size_t n = x.size();

    // sort x and weights in x order
    auto perm = utils::get_order(x);
    auto xx = x;
    auto w = weights;
    for (size_t i = 0; i < n; i++) {
        xx[i] = x[perm[i]];
        if (w.size() > 0)
            w[i] = weights[perm[i]];
    }

    // compute weighted ranks and the "average rank" (corresponds to the
    // median)
    auto ranks = rank(xx, w, "average");
    if (weights.size() == 0)
        weights = std::vector<double>(n, 1.0);
    double rank_avrg = utils::perm_sum(weights, 2) / utils::sum(weights);

    // weighted median splits data below and above rank_avrg
    size_t i = 0;
    while (ranks[i] < rank_avrg)
        i++;
    if (ranks[i] == rank_avrg)
        return xx[i];
    else
        return 0.5 * (xx[i - 1] + xx[i]);
}
}
}
