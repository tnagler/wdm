#include <cmath>
#include <random>
#include <vector>
#include <wdm.hpp>

#include "test_helpers.hpp"

namespace {

using test::check;
using test::check_near;
using test::check_throws;
using wdm::impl::default_tie_seeds;

double
continuous_xi_for_inference(const std::vector<double>& y,
                            std::vector<double> probabilities)
{
  double probability_sum = wdm::utils::sum(probabilities);
  for (auto& probability : probabilities)
    probability /= probability_sum;
  auto weighted_ranks = wdm::impl::rank0(y, probabilities, "max");
  double weighted_edge_difference = 0.0;
  for (size_t i = 0; i + 1 < weighted_ranks.size(); ++i) {
    weighted_edge_difference +=
      probabilities[i] * std::fabs(weighted_ranks[i + 1] - weighted_ranks[i]);
  }
  return 1.0 - 3.0 * weighted_edge_difference;
}

void
check_cxi_null_variance(const std::vector<double>& weights,
                        const std::string& what)
{
  std::vector<double> x(weights.size()), y(weights.size());
  for (size_t i = 0; i < weights.size(); ++i) {
    x[i] = static_cast<double>(i);
    y[i] = static_cast<double>(i);
  }

  auto inference = wdm::impl::cxi(x, y, weights, true);
  std::mt19937 generator(20260821);
  double squared_deviation_sum = 0.0;
  size_t replications = 5000;
  for (size_t replication = 0; replication < replications; ++replication) {
    std::shuffle(y.begin(), y.end(), generator);
    double deviation =
      continuous_xi_for_inference(y, weights) - std::get<2>(inference);
    squared_deviation_sum += deviation * deviation;
  }

  double empirical_variance = squared_deviation_sum / replications;
  double analytic_variance = std::get<1>(inference) * std::get<1>(inference);
  check(std::fabs(empirical_variance / analytic_variance - 1.0) < 0.12, what);
}

//! for strictly increasing data the ranks are 1, ..., n, so the numerator is
//! n - 1 and xi is 1 - 3 / (n + 1).
void
test_cxi()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6, 7, 8 };
  double n = static_cast<double>(x.size());
  check_near(
    wdm::wdm(x, x, "cxi"), 1.0 - 3.0 / (n + 1.0), "xi of increasing data");
  check_near(wdm::wdm(x, x, "chatterjee"),
             wdm::wdm(x, x, "cxi"),
             "xi method names are aliases");
  check_near(wdm::wdm(x, x, "xi"),
             wdm::wdm(x, x, "cxi"),
             "short xi method name is an alias");

  // xi picks up a non-monotonic relationship, and is asymmetric: the square is
  // a measurable function of the argument, but not the other way around
  std::vector<double> v, v_sq;
  for (int i = -20; i <= 20; ++i) {
    if (i == 0)
      continue;
    v.push_back(i);
    v_sq.push_back(i * i);
  }
  check(wdm::wdm(v, v_sq, "cxi") > 0.8, "xi detects a non-monotonic function");
  check(wdm::wdm(v, v_sq, "cxi") > wdm::wdm(v_sq, v, "cxi") + 0.5,
        "xi is asymmetric");

  // Reordering observations does not matter when the predictor has no ties.
  std::vector<double> v_rev(v.rbegin(), v.rend());
  std::vector<double> v_sq_rev(v_sq.rbegin(), v_sq.rend());
  check_near(wdm::wdm(v_rev, v_sq_rev, "cxi"),
             wdm::wdm(v, v_sq, "cxi"),
             "xi ignores the order of the observations");

  // uniform weights must not change anything
  check_near(wdm::wdm(v, v_sq, "cxi", std::vector<double>(v.size(), 1.0)),
             wdm::wdm(v, v_sq, "cxi"),
             "xi with uniform weights");

  // Unequal weights use the base point of each edge. For this example the
  // weighted numerator is 1/4 and the weighted denominator is 31/108.
  std::vector<double> short_x{ 1, 2, 3 };
  std::vector<double> short_y{ 1, 3, 2 };
  std::vector<double> unequal_weights{ 1, 2, 3 };
  check_near(wdm::wdm(short_x, short_y, "cxi", unequal_weights),
             4.0 / 31.0,
             "weighted xi uses base-point edge weights");
  for (auto& weight : unequal_weights)
    weight *= 10.0;
  check_near(wdm::wdm(short_x, short_y, "cxi", unequal_weights),
             4.0 / 31.0,
             "weighted xi is invariant to weight scaling");

  auto unequal_inference = wdm::impl::cxi(short_x, short_y, { 1, 2, 3 }, true);
  auto scaled_inference =
    wdm::impl::cxi(short_x, short_y, { 10, 20, 30 }, true);
  check_near(std::get<1>(unequal_inference),
             std::sqrt(7.0 / 120.0),
             "weighted xi uses the full conditional null variance");
  check_near(std::get<2>(unequal_inference),
             23.0 / 72.0,
             "weighted xi has the finite-sample null mean");
  check_near(std::get<3>(unequal_inference),
             1.0 / 4.0,
             "weighted xi exposes the continuous inferential estimate");
  check_near(std::get<1>(scaled_inference),
             std::get<1>(unequal_inference),
             "weighted xi standard error is invariant to weight scaling");
  check_near(std::get<2>(scaled_inference),
             std::get<2>(unequal_inference),
             "weighted xi null mean is invariant to weight scaling");
  check_near(std::get<3>(scaled_inference),
             std::get<3>(unequal_inference),
             "weighted xi inferential estimate is invariant to weight scaling");

  std::vector<double> large_x(1000), large_y(1000);
  for (size_t i = 0; i < large_x.size(); ++i) {
    large_x[i] = static_cast<double>(i);
    large_y[i] = static_cast<double>((37 * i) % large_x.size());
  }
  auto equal_inference = wdm::impl::cxi(large_x, large_y, {}, true);
  check(std::fabs(std::get<1>(equal_inference) /
                    std::sqrt(2.0 / (5.0 * large_x.size())) -
                  1.0) < 0.002,
        "equal-weight xi standard error approaches sqrt(2 / (5 n))");
  check_near(std::get<2>(equal_inference),
             1.0 / (large_x.size() * large_x.size()),
             "equal-weight xi has the finite-sample null mean");

  std::vector<double> smooth_weights(200), alternating_weights(200);
  for (size_t i = 0; i < smooth_weights.size(); ++i) {
    smooth_weights[i] = 1.0 + static_cast<double>(i) / smooth_weights.size();
    alternating_weights[i] = (i % 2 == 0) ? 0.5 : 1.5;
  }
  check_cxi_null_variance(
    smooth_weights,
    "full xi null variance agrees with simulation for smooth weights");
  check_cxi_null_variance(
    alternating_weights,
    "full xi null variance agrees with simulation for alternating weights");

  check_near(wdm::wdm({ 1, 2, 3, 4 }, { 1, 2, 2, 1 }, "cxi"),
             0.0,
             "xi uses the general denominator for tied responses");

  check_throws(
    [&]() {
      wdm::wdm(short_x, short_y, "cxi", { 1, -1, 1 });
    },
    "xi rejects negative weights");
  check_throws(
    [&]() {
      wdm::wdm(short_x, short_y, "cxi", { 1, INFINITY, 1 });
    },
    "xi rejects nonfinite weights");
  check_throws(
    [&]() {
      wdm::wdm(short_x, short_y, "cxi", { 0, 0, 0 });
    },
    "xi rejects zero total weight");
  check_throws(
    [&]() {
      wdm::wdm(short_x, { 1, 1, 1 }, "cxi");
    },
    "xi rejects a constant response");

  std::vector<double> tied_predictor{ 2, 1, 1, 2, 1 };
  std::vector<double> tied_predictor_response{ 10, 20, 30, 40, 50 };
  std::vector<double> tied_predictor_weights{ 1, 2, 3, 4, 5 };
  std::vector<int> tie_seeds{ 17, 29, 43 };
  check_near(wdm::wdm(tied_predictor,
                      tied_predictor_response,
                      "cxi",
                      tied_predictor_weights,
                      true,
                      tie_seeds),
             wdm::wdm(tied_predictor,
                      tied_predictor_response,
                      "cxi",
                      tied_predictor_weights,
                      true,
                      tie_seeds),
             "seeded predictor-tie breaking is reproducible");
  check_near(wdm::wdm(v, v_sq, "cxi", {}, true, { 1 }),
             wdm::wdm(v, v_sq, "cxi", {}, true, { 2 }),
             "predictor-tie seeds do not affect untied data");

  // Without a seed the tie ordering must still be reproducible: two large tie
  // groups leave enough orderings that a random seed would essentially never
  // repeat itself.
  std::vector<double> wide_x(60), wide_y(60);
  for (size_t i = 0; i < wide_x.size(); ++i) {
    wide_x[i] = static_cast<double>(i / 30);
    wide_y[i] = static_cast<double>((37 * i) % wide_x.size());
  }
  double unseeded = wdm::wdm(wide_x, wide_y, "cxi");
  for (int repetition = 0; repetition < 5; ++repetition) {
    check_near(wdm::wdm(wide_x, wide_y, "cxi"),
               unseeded,
               "unseeded predictor-tie breaking is reproducible");
  }
  check_near(wdm::wdm(wide_x, wide_y, "cxi", {}, true, default_tie_seeds()),
             unseeded,
             "the unseeded default is default_tie_seeds()");
  wdm::Indep_test unseeded_test(wide_x, wide_y, "cxi");
  wdm::Indep_test unseeded_test_again(wide_x, wide_y, "cxi");
  check_near(unseeded_test.p_value(),
             unseeded_test_again.p_value(),
             "unseeded xi inference is reproducible");

  // The randomization is still available: explicit seeds must reach it.
  std::vector<double> seeded_values;
  for (int seed = 1; seed <= 6; ++seed)
    seeded_values.push_back(
      wdm::wdm(wide_x, wide_y, "cxi", {}, true, { seed }));
  bool some_seed_differs = false;
  for (const auto& value : seeded_values)
    some_seed_differs = some_seed_differs || value != seeded_values.front();
  check(some_seed_differs, "explicit seeds change the tie ordering");

  std::vector<double> shifted_response = tied_predictor_response;
  for (auto& response : shifted_response)
    response += 100.0;
  std::vector<double> sorted_predictor = tied_predictor;
  std::vector<double> sorted_response = tied_predictor_response;
  std::vector<double> sorted_weights = tied_predictor_weights;
  wdm::impl::sort_chatterjee_observations(
    sorted_predictor, sorted_response, sorted_weights, tie_seeds);
  wdm::impl::sort_chatterjee_observations(
    tied_predictor, shifted_response, tied_predictor_weights, tie_seeds);
  check(std::is_sorted(sorted_predictor.begin(), sorted_predictor.end()),
        "Chatterjee observations are sorted by the predictor");
  for (size_t i = 0; i < sorted_response.size(); ++i) {
    check_near(sorted_response[i],
               10.0 * sorted_weights[i],
               "predictor-tie breaking keeps weights with observations");
    check_near(shifted_response[i],
               sorted_response[i] + 100.0,
               "predictor-tie breaking does not consult the response");
  }

  auto tied_predictor_inference = wdm::impl::cxi(
    sorted_predictor, sorted_response, sorted_weights, true, "max", tie_seeds);
  wdm::Indep_test tied_predictor_test(sorted_predictor,
                                      sorted_response,
                                      "cxi",
                                      sorted_weights,
                                      true,
                                      "two-sided",
                                      tie_seeds);
  check_near(tied_predictor_test.estimate(),
             std::get<0>(tied_predictor_inference),
             "xi test reuses the seeded predictor-tie ordering");
  check_near(tied_predictor_test.statistic(),
             (std::get<3>(tied_predictor_inference) -
              std::get<2>(tied_predictor_inference)) /
               std::get<1>(tied_predictor_inference),
             "xi inference uses the realized predictor-tie ordering");

  wdm::Indep_test test(v, v_sq, "cxi");
  check(std::isfinite(test.p_value()), "xi p-value is finite");

  wdm::Indep_test weighted_test(short_x, short_y, "cxi", { 1, 2, 3 });
  check_near(weighted_test.estimate(),
             std::get<0>(unequal_inference),
             "xi test reports the denominator-corrected estimate");
  check_near(weighted_test.statistic(),
             (std::get<3>(unequal_inference) - std::get<2>(unequal_inference)) /
               std::get<1>(unequal_inference),
             "xi test uses the derived continuous-response statistic");
  wdm::Indep_test scaled_weight_test(short_x, short_y, "cxi", { 10, 20, 30 });
  check_near(scaled_weight_test.statistic(),
             weighted_test.statistic(),
             "xi test statistic is invariant to weight scaling");
  check_near(scaled_weight_test.p_value(),
             weighted_test.p_value(),
             "xi p-value is invariant to weight scaling");

  auto unweighted_discrete_inference =
    wdm::impl::cxi(short_x, short_y, {}, true, "max", {}, false);
  wdm::Indep_test unweighted_discrete_test(
    short_x, short_y, "cxi", {}, true, "two-sided", {}, false);
  check(std::isfinite(unweighted_discrete_test.p_value()),
        "unweighted discrete-response xi inference remains available");
  check_near(unweighted_discrete_test.statistic(),
             (std::get<3>(unweighted_discrete_inference) -
              std::get<2>(unweighted_discrete_inference)) /
               std::get<1>(unweighted_discrete_inference),
             "declared discrete response uses tied-response inference");
  check_throws(
    [&]() {
      wdm::Indep_test weighted_discrete_test(
        short_x, short_y, "cxi", { 1, 2, 3 }, true, "two-sided", {}, false);
    },
    "unequally weighted discrete-response xi inference is unavailable");

  std::vector<double> tied_x{ 1, 2, 3, 4 };
  std::vector<double> tied_y{ 1, 2, 2, 1 };
  wdm::Indep_test tied_test(tied_x, tied_y, "cxi");
  check(std::isfinite(tied_test.p_value()),
        "unweighted tied-response xi inference remains available");
  wdm::Indep_test uniformly_weighted_tied_test(
    tied_x, tied_y, "cxi", { 10, 10, 10, 10 });
  check_near(uniformly_weighted_tied_test.statistic(),
             tied_test.statistic(),
             "uniform weights retain tied-response xi inference");
  check_throws(
    [&]() {
      wdm::Indep_test weighted_tied_test(tied_x, tied_y, "cxi", { 1, 2, 1, 2 });
    },
    "weighted tied-response xi inference is unavailable");
}

} // namespace

int
main()
{
  test_cxi();
  return test::finish();
}
