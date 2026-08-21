// Copyright © 2020 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdm/blob/master/LICENSE.

#pragma once

#include "wdm/bbeta.hpp"
#include "wdm/cxi.hpp"
#include "wdm/hoeffd.hpp"
#include "wdm/ktau.hpp"
#include "wdm/methods.hpp"
#include "wdm/nan_handling.hpp"
#include "wdm/prho.hpp"
#include "wdm/srho.hpp"

//! Weighted dependence measures
namespace wdm {

//! Calculates a weighted or unweighted dependence measure.
//! @param x, y input vectors of equal length. For Chatterjee's xi, `x` is the
//!   predictor and `y` is the response.
//! @param method the dependence measure; see details for possible values.
//! @param weights optional case weights. Nonempty weights must match the input
//!   length, be finite and nonnegative, and have a positive sum. Their scale
//!   does not affect the result, and zero-weight rows are ignored.
//! @param remove_missing if `true`, rows containing a `NaN` are removed;
//!   otherwise a `NaN` raises an exception.
//! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
//!
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$
//!   - `"chatterjee"`, `"cxi"`, `"xi"`: Chatterjee's \f$ \xi \f$
//!
//! @note Chatterjee's xi is asymmetric. Predictor ties are broken without
//!   consulting the response; use `seeds` for a reproducible ordering.
//! @return The requested estimate, or `NaN` if missing-value removal leaves
//!   fewer than two observations (five for Hoeffding's D).
//! @throws std::runtime_error for size mismatches, invalid weights, unknown
//!   methods, disallowed missing or insufficient input, or an undefined
//!   Chatterjee estimate with a constant response.
inline double
wdm(std::vector<double> x,
    std::vector<double> y,
    std::string method,
    std::vector<double> weights = std::vector<double>(),
    bool remove_missing = true,
    std::vector<int> seeds = std::vector<int>())
{
  utils::check_sizes(x, y, weights);
  // na handling
  if (utils::preproc(x, y, weights, method, remove_missing) == "return_nan")
    return std::numeric_limits<double>::quiet_NaN();

  if (methods::is_hoeffding(method))
    return impl::hoeffd(x, y, weights);
  if (methods::is_kendall(method))
    return impl::ktau(x, y, weights);
  if (methods::is_pearson(method))
    return impl::prho(x, y, weights);
  if (methods::is_spearman(method))
    return impl::srho(x, y, weights);
  if (methods::is_blomqvist(method))
    return impl::bbeta(x, y, weights);
  if (methods::is_chatterjee(method)) {
    auto xi_and_std = impl::cxi(x, y, weights, false, "max", seeds);
    return std::get<0>(xi_and_std);
  }
  throw std::runtime_error("method not implemented.");
}

//! Asymptotic independence test based on a dependence measure.
//!
//! The test stores the estimate, transformed test statistic, effective sample
//! size, and p-value. Weighted transformations use Kish's effective sample
//! size. The approximation must have enough effective observations for the
//! selected method.
//!
//! @details
//! Available methods:
//!   - `"pearson"`, `"prho"`, `"cor"`: Pearson correlation
//!   - `"spearman"`, `"srho"`, `"rho"`: Spearman's \f$ \rho \f$
//!   - `"kendall"`, `"ktau"`, `"tau"`: Kendall's \f$ \tau \f$
//!   - `"blomqvist"`, `"bbeta"`, `"beta"`: Blomqvist's \f$ \beta \f$
//!   - `"hoeffding"`, `"hoeffd"`, `"d"`: Hoeffding's \f$ D \f$
//!   - `"chatterjee"`, `"cxi"`, `"xi"`: Chatterjee's \f$ \xi \f$
//!
//! Hoeffding's D supports only the two-sided alternative. Other methods support
//! `"two-sided"`, `"less"`, and `"greater"`.
//!
//! @note Weighted analytic inference for Chatterjee's xi assumes that the
//!   weights are fixed or depend only on `x`, the normalized weights are
//!   diffuse, and `y` is continuous. The weighted estimate remains available
//!   when `y` is tied, but analytic inference with unequal weights does not.
//!   For a continuous response, `estimate()` reports the general
//!   denominator-corrected coefficient while `statistic()` standardizes the
//!   analytically covered approximation \f$1 - 3 A\f$.
//!
class Indep_test
{
public:
  Indep_test() = delete;

  //! Constructs and evaluates an independence test.
  //! @param x, y input vectors of equal length.
  //! @param method the dependence measure; see class details for possible
  //! values.
  //! @param weights optional finite, nonnegative case weights with positive
  //!   total weight.
  //! @param remove_missing if `true`, rows containing a `NaN` are removed;
  //!   otherwise a `NaN` raises an exception.
  //! @param alternative indicates the alternative hypothesis and must be one
  //!    of `"two-sided"`, `"greater"` or `"less"`; `"greater"` corresponds
  //!    to positive association, `"less"` to negative association. For
  //!    Hoeffding's \f$ D \f$, only `"two-sided"` is allowed. The natural
  //!    one-sided alternative for Chatterjee's xi is `"greater"`.
  //! @param seeds optional seeds for random Chatterjee predictor-tie breaking.
  //! @param y_continuous whether the Chatterjee response distribution is known
  //!    to be continuous. Set this to `false` for a discrete response even if
  //!    the sample has no observed response ties. Observed ties always override
  //!    this value.
  //! @throws std::runtime_error for invalid inputs, method or alternative
  //!   names, unsupported Hoeffding alternatives, or unavailable weighted
  //!   Chatterjee inference for a discrete or tied response.
  Indep_test(std::vector<double> x,
             std::vector<double> y,
             std::string method,
             std::vector<double> weights = std::vector<double>(),
             bool remove_missing = true,
             std::string alternative = "two-sided",
             std::vector<int> seeds = std::vector<int>(),
             bool y_continuous = true)
    : method_(method)
    , alternative_(alternative)
  {
    utils::check_sizes(x, y, weights);
    if (utils::preproc(x, y, weights, method, remove_missing) == "return_nan") {
      n_eff_ = utils::effective_sample_size(x.size(), weights);
      estimate_ = std::numeric_limits<double>::quiet_NaN();
      statistic_ = std::numeric_limits<double>::quiet_NaN();
      p_value_ = std::numeric_limits<double>::quiet_NaN();
    } else {
      n_eff_ = utils::effective_sample_size(x.size(), weights);
      if (methods::is_chatterjee(method)) {
        auto stats = impl::cxi(x, y, weights, true, "max", seeds, y_continuous);
        estimate_ = std::get<0>(stats);
        statistic_ =
          (std::get<3>(stats) - std::get<2>(stats)) / std::get<1>(stats);
      } else {
        estimate_ = wdm(x, y, method, weights, false);
        statistic_ =
          compute_test_stat(estimate_, method, n_eff_, x, y, weights);
      }
      p_value_ = compute_p_value(statistic_, method, alternative, n_eff_);
    }
  }

  //! Returns the requested method name.
  std::string method() const { return method_; }

  //! Returns the requested alternative hypothesis.
  std::string alternative() const { return alternative_; }

  //! Returns Kish's effective sample size after missing-value removal.
  double n_eff() const { return n_eff_; }

  //! Returns the estimated dependence measure.
  double estimate() const { return estimate_; }

  //! Returns the method-specific transformed test statistic.
  double statistic() const { return statistic_; }

  //! Returns the asymptotic p-value.
  double p_value() const { return p_value_; }

private:
  inline double compute_test_stat(double estimate,
                                  std::string method,
                                  double n_eff,
                                  const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& weights)
  {
    // prevent overflow in atanh
    if (estimate >= 1.0)
      estimate = 1 - 1e-12;
    if (estimate <= -1.0)
      estimate = -1 + 1e-12;

    double stat;
    if (methods::is_hoeffding(method)) {
      stat = estimate / 30.0 + 1.0 / (36.0 * n_eff);
    } else if (methods::is_kendall(method)) {
      stat = estimate * impl::ktau_stat_adjust(x, y, weights);
    } else if (methods::is_pearson(method)) {
      stat = std::atanh(estimate) * std::sqrt(n_eff - 3);
    } else if (methods::is_spearman(method)) {
      stat = std::atanh(estimate) * std::sqrt((n_eff - 3) / 1.06);
    } else if (methods::is_blomqvist(method)) {
      stat = std::atanh(estimate) * std::sqrt(n_eff);
    } else {
      throw std::runtime_error("method not implemented.");
    }

    return stat;
  }

  inline double compute_p_value(double statistic,
                                std::string method,
                                std::string alternative,
                                double n_eff = 0.0)
  {
    double p_value;
    if (methods::is_hoeffding(method)) {
      if (n_eff == 0.0)
        throw std::runtime_error("must provide n_eff for method 'hoeffd'.");
      if (alternative != "two-sided")
        throw std::runtime_error(
          "only two-sided test available for Hoeffding's D.");
      p_value = impl::phoeffb(statistic, n_eff);
    } else {
      if (alternative == "two-sided") {
        p_value = 2 * utils::normalCDF(-std::abs(statistic));
      } else if (alternative == "less") {
        p_value = utils::normalCDF(statistic);
      } else if (alternative == "greater") {
        p_value = 1 - utils::normalCDF(statistic);
      } else {
        throw std::runtime_error("alternative not implemented.");
      }
    }

    return p_value;
  }

  std::string method_;
  std::string alternative_;
  double n_eff_;
  double estimate_;
  double statistic_;
  double p_value_;
};

}
