/*
int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/

// #include <wdm/include/wdm.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <wdm.hpp>

namespace {

int failures = 0;

void
check(bool condition, const std::string& what)
{
  if (!condition) {
    std::cout << "FAILED: " << what << std::endl;
    ++failures;
  }
}

void
check_near(double actual, double expected, const std::string& what)
{
  if (!(std::fabs(actual - expected) < 1e-10)) {
    std::cout << "FAILED: " << what << " (got " << std::setprecision(12)
              << actual << ", expected " << expected << ")" << std::endl;
    ++failures;
  }
}

template<typename Function>
void
check_throws(Function function, const std::string& what)
{
  bool threw = false;
  try {
    function();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, what);
}

bool
all_close(std::vector<double> x, std::vector<double> y, double tol = 1e-12)
{
  if (x.size() != y.size())
    return false;
  for (size_t i = 0; i < x.size(); ++i) {
    if (std::fabs(x[i] - y[i]) > tol)
      return false;
  }
  return true;
}

std::vector<double>
sorted(std::vector<double> x)
{
  std::sort(x.begin(), x.end());
  return x;
}

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

//! successive differences, the first one taken from 0.
std::vector<double>
diff(const std::vector<double>& x)
{
  std::vector<double> res(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    res[i] = x[i] - (i > 0 ? x[i - 1] : 0.0);
  return res;
}

//! every tie breaking method assigns the same ranks to a tie group, differing
//! only in how it distributes them; `"random"` draws the order.
void
test_rank_ties()
{
  std::vector<int> seeds{ 1, 2, 3 };

  // a tie group of six, occupying the rank slots 2, ..., 7
  std::vector<double> x{ 1, 2, 2, 2, 2, 2, 2, 3, 4, 5 };
  std::vector<double> consecutive{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

  auto first = wdm::impl::rank(x, {}, "first");
  check(all_close(first, consecutive),
        "'first' ranks tied values in order of appearance");
  check(
    all_close(wdm::impl::rank(x, {}, "min"), { 1, 2, 2, 2, 2, 2, 2, 8, 9, 10 }),
    "'min' assigns tied values the minimum rank");
  check(all_close(wdm::impl::rank(x, {}, "average"),
                  { 1, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 8, 9, 10 }),
        "'average' assigns tied values the average rank");

  // 'random' assigns the same ranks as 'first', in a different order
  auto random = wdm::impl::rank(x, {}, "random", seeds);
  check(all_close(sorted(random), consecutive),
        "'random' permutes the ranks assigned by 'first'");
  check(all_close(random, wdm::impl::rank(x, {}, "random", seeds)),
        "'random' is reproducible for identical seeds");

  // weights are normalized to mean one, so a single tie group spanning the
  // whole vector receives the cumulative weights in the order it drew
  std::vector<double> tied(5, 2.0);
  std::vector<double> weights{ 1, 2, 3, 4, 5 };
  auto weighted = sorted(wdm::impl::rank(tied, weights, "random", seeds));
  check(
    all_close(sorted(diff(weighted)), { 1. / 3, 2. / 3, 1., 4. / 3, 5. / 3 }),
    "weighted 'random' ranks are the cumulative weights in drawn order");
  check_near(weighted.back(), 5.0, "largest rank of a tie group");

  // separate tie groups are shuffled independently
  std::vector<double> pairs;
  for (size_t i = 0; i < 20; ++i)
    pairs.insert(pairs.end(), 2, static_cast<double>(i));
  auto shuffled = wdm::impl::rank(pairs, {}, "random", seeds);
  bool ascending = false, descending = false;
  for (size_t i = 0; i < shuffled.size(); i += 2) {
    ascending = ascending || (shuffled[i] < shuffled[i + 1]);
    descending = descending || (shuffled[i] > shuffled[i + 1]);
  }
  check(ascending && descending, "tie groups are shuffled independently");
  check(all_close(sorted(diff(sorted(shuffled))),
                  std::vector<double>(pairs.size(), 1.0)),
        "'random' leaves no ties behind");

  // NaNs keep their place and take no rank
  std::vector<double> with_nan{ 2, NAN, 2, 2 };
  auto nan_ranks = wdm::impl::rank(with_nan, {}, "random", seeds);
  check(std::isnan(nan_ranks[1]), "'random' preserves NaNs");
  nan_ranks.erase(nan_ranks.begin() + 1);
  check(all_close(sorted(nan_ranks), { 1, 2, 3 }),
        "'random' ranks the non-NaN observations among themselves");
}

//! a rank under `"min"` is the weight of the strictly smaller observations, one
//! under `"max"` the weight of those less than or equal, so the two differ by
//! the weight of the tied batch.
void
test_rank0_ties_methods()
{
  std::vector<double> x{ 1, 3, 2, 5, 3, 2, 20, 15 };
  std::vector<double> expected_min{ 0, 3, 1, 5, 3, 1, 7, 6 };
  std::vector<double> expected_max{ 1, 5, 3, 6, 5, 3, 8, 7 };

  auto r_min = wdm::impl::rank0(x);
  auto r_max = wdm::impl::rank0(x, std::vector<double>(), "max");
  for (size_t i = 0; i < x.size(); ++i) {
    check_near(r_min[i], expected_min[i], "rank0 'min'");
    check_near(r_max[i], expected_max[i], "rank0 'max'");
  }

  bool threw = false;
  try {
    wdm::impl::rank0(x, std::vector<double>(), "foo");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "rank0 rejects an unknown ties_method");
}

//! rank0 takes the weights as they are, so an average rank sits halfway up the
//! weight of its own tie group.
void
test_rank0_average_ties()
{
  std::vector<double> x{ 1, 2, 2, 2, 2, 2, 2, 3, 4, 5 };
  std::vector<double> weights{ 1, 1, 2, 2, 1, 3, 1, 1, 1, 1 };

  check(all_close(wdm::impl::rank0(x, weights),
                  { 0, 1, 1, 1, 1, 1, 1, 11, 12, 13 }),
        "weighted rank0 'min' accumulates the weights");
  check(all_close(wdm::impl::rank0(x, {}, "average"),
                  { 0, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 7, 8, 9 }),
        "rank0 'average' assigns tied values the average rank");
  check(all_close(wdm::impl::rank0(x, weights, "average"),
                  { 0, 5, 5, 5, 5, 5, 5, 11, 12, 13 }),
        "weighted rank0 'average' assigns the average weighted rank");
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

  // reordering the observations must not change either direction
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
  check_near(std::get<1>(scaled_inference),
             std::get<1>(unequal_inference),
             "weighted xi standard error is invariant to weight scaling");
  check_near(std::get<2>(scaled_inference),
             std::get<2>(unequal_inference),
             "weighted xi null mean is invariant to weight scaling");

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

  check_throws([&]() { wdm::wdm(short_x, short_y, "cxi", { 1, -1, 1 }); },
               "xi rejects negative weights");
  check_throws([&]() { wdm::wdm(short_x, short_y, "cxi", { 1, INFINITY, 1 }); },
               "xi rejects nonfinite weights");
  check_throws([&]() { wdm::wdm(short_x, short_y, "cxi", { 0, 0, 0 }); },
               "xi rejects zero total weight");
  check_throws([&]() { wdm::wdm(short_x, { 1, 1, 1 }, "cxi"); },
               "xi rejects a constant response");

  wdm::Indep_test test(v, v_sq, "cxi");
  check(std::isfinite(test.p_value()), "xi p-value is finite");

  wdm::Indep_test weighted_test(short_x, short_y, "cxi", { 1, 2, 3 });
  check_near(weighted_test.statistic(),
             (std::get<0>(unequal_inference) - std::get<2>(unequal_inference)) /
               std::get<1>(unequal_inference),
             "xi test statistic uses the finite-sample null mean");
  wdm::Indep_test scaled_weight_test(short_x, short_y, "cxi", { 10, 20, 30 });
  check_near(scaled_weight_test.statistic(),
             weighted_test.statistic(),
             "xi test statistic is invariant to weight scaling");
  check_near(scaled_weight_test.p_value(),
             weighted_test.p_value(),
             "xi p-value is invariant to weight scaling");

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
  // input vectors
  std::vector<double> x{ 1, 3, 2, 5, 3, 2, 20, 15 };
  std::vector<double> y{ 2, 12, 4, 7, 8, 14, 17, 6 };

  // weights
  std::vector<double> w{ 1, 1, 2, 2, 1, 0, 0.5, 0.3 };

  std::cout << "unweighted Kendall's tau: " << wdm::wdm(x, y, "kendall")
            << std::endl;
  std::cout << "weighted Kendall's tau: " << wdm::wdm(x, y, "kendall", w)
            << std::endl;

  // weighted independence test
  wdm::Indep_test test(x, y, "kendall", w);
  std::cout << "statistic: " << test.statistic() << std::endl;
  std::cout << "p-value: " << test.p_value() << std::endl;

  test_rank_ties();
  test_rank0_ties_methods();
  test_rank0_average_ties();
  test_cxi();

  if (failures > 0) {
    std::cout << failures << " check(s) failed" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}
