#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <wdm.hpp>

#include "test_helpers.hpp"

namespace {

void
test_statistic_transformations()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<double> y{ 2, 6, 4, 3, 7, 1, 8, 5 };

  wdm::Indep_test pearson(x, y, "pearson");
  test::check_near(pearson.statistic(),
                   std::atanh(pearson.estimate()) * std::sqrt(5.0),
                   "Pearson Fisher transformation");

  wdm::Indep_test spearman(x, y, "spearman");
  test::check_near(spearman.statistic(),
                   std::atanh(spearman.estimate()) * std::sqrt(5.0 / 1.06),
                   "Spearman Fisher transformation");

  wdm::Indep_test kendall(x, y, "kendall");
  test::check_near(kendall.statistic(),
                   kendall.estimate() * wdm::impl::ktau_stat_adjust(x, y, {}),
                   "Kendall tie-adjusted transformation");

  wdm::Indep_test blomqvist(x, y, "blomqvist");
  test::check_near(blomqvist.statistic(),
                   std::atanh(blomqvist.estimate()) * std::sqrt(8.0),
                   "Blomqvist Fisher transformation");

  wdm::Indep_test hoeffding(x, y, "hoeffding");
  test::check_near(hoeffding.statistic(),
                   hoeffding.estimate() / 30.0 + 1.0 / (36.0 * 8.0),
                   "Hoeffding B transformation");

  std::vector<double> decreasing(x.rbegin(), x.rend());
  wdm::Indep_test negative_endpoint(x, decreasing, "pearson");
  test::check(negative_endpoint.statistic() < -10.0,
              "negative Pearson endpoint remains negative after clamping");
  test::check(negative_endpoint.p_value() < 1e-10,
              "negative Pearson endpoint has a small two-sided p-value");
}

void
test_hoeffding_tail_table()
{
  test::check_near(wdm::utils::linear_interp(1.25, { 1, 2 }, { 10, 20 }),
                   12.5,
                   "linear interpolation uses the lower endpoint weight");
  test::check_near(wdm::impl::phoeffb(2.2 / std::pow(wdm::impl::pi, 4), 2.0),
                   0.5297,
                   "Hoeffding table lower boundary");
  test::check_near(wdm::impl::phoeffb(11.0 / std::pow(wdm::impl::pi, 4), 2.0),
                   0.0025,
                   "Hoeffding table value at 5.5");
  test::check_near(wdm::impl::phoeffb(11.5 / std::pow(wdm::impl::pi, 4), 2.0),
                   0.00195,
                   "Hoeffding interpolation between table nodes");
  test::check_near(wdm::impl::phoeffb(12.0 / std::pow(wdm::impl::pi, 4), 2.0),
                   0.0014,
                   "Hoeffding table value at 6.0");
  test::check_near(wdm::impl::phoeffb(17.0 / std::pow(wdm::impl::pi, 4), 2.0),
                   0.0001,
                   "Hoeffding table upper boundary");
}

void
test_alternatives_and_p_values()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<double> y{ 2, 1, 4, 3, 7, 5, 8, 6 };
  for (const auto& method : std::vector<std::string>{
         "pearson", "spearman", "kendall", "blomqvist", "chatterjee" }) {
    wdm::Indep_test two_sided(x, y, method);
    wdm::Indep_test less(x, y, method, {}, true, "less");
    wdm::Indep_test greater(x, y, method, {}, true, "greater");
    test::check_near(less.p_value() + greater.p_value(),
                     1.0,
                     method + " one-sided p-values are complementary");
    test::check_near(two_sided.p_value(),
                     2.0 * std::min(less.p_value(), greater.p_value()),
                     method + " two-sided p-value doubles the smaller tail");
    test::check_near(greater.p_value(),
                     1.0 - wdm::utils::normalCDF(two_sided.statistic()),
                     method + " greater p-value uses the upper normal tail");
    test::check(two_sided.alternative() == "two-sided",
                method + " reports its alternative");
  }

  test::check_throws(
    [&]() { wdm::Indep_test(x, y, "pearson", {}, true, "unknown"); },
    "inference rejects an unknown alternative");
  test::check_throws(
    [&]() { wdm::Indep_test(x, y, "hoeffding", {}, true, "greater"); },
    "Hoeffding inference rejects one-sided alternatives");
}

void
test_aliases()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<double> y{ 2, 1, 4, 3, 7, 5, 8, 6 };
  std::vector<std::vector<std::string>> aliases{
    { "pearson", "prho", "cor" },   { "spearman", "srho", "rho" },
    { "kendall", "ktau", "tau" },   { "blomqvist", "bbeta", "beta" },
    { "hoeffding", "hoeffd", "d" }, { "chatterjee", "cxi", "xi" }
  };
  for (const auto& method_aliases : aliases) {
    wdm::Indep_test expected(x, y, method_aliases.front());
    for (const auto& method : method_aliases) {
      wdm::Indep_test result(x, y, method);
      test::check_near(result.estimate(),
                       expected.estimate(),
                       method + " inference estimate alias");
      test::check_near(result.statistic(),
                       expected.statistic(),
                       method + " inference statistic alias");
      test::check_near(result.p_value(),
                       expected.p_value(),
                       method + " inference p-value alias");
      test::check(result.method() == method,
                  method + " preserves the requested method name");
    }
  }
}

double
kendall_score(const std::vector<double>& x, const std::vector<double>& y)
{
  double score = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = i + 1; j < x.size(); ++j) {
      if (x[i] != x[j] && y[i] != y[j])
        score += (x[i] < x[j]) == (y[i] < y[j]) ? 1.0 : -1.0;
    }
  }
  return score;
}

std::vector<size_t>
tie_group_sizes(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  std::vector<size_t> group_sizes;
  for (size_t begin = 0, end; begin < values.size(); begin = end) {
    end = begin + 1;
    while (end < values.size() && values[end] == values[begin])
      ++end;
    if (end - begin > 1)
      group_sizes.push_back(end - begin);
  }
  return group_sizes;
}

double
kendall_score_variance(const std::vector<double>& x,
                       const std::vector<double>& y)
{
  double sample_size = static_cast<double>(x.size());
  double x_pair_term = 0.0;
  double y_pair_term = 0.0;
  double x_triplet_term = 0.0;
  double y_triplet_term = 0.0;
  double x_variance_term = 0.0;
  double y_variance_term = 0.0;
  for (const auto& group_size : tie_group_sizes(x)) {
    double size = static_cast<double>(group_size);
    x_pair_term += size * (size - 1.0);
    x_triplet_term += size * (size - 1.0) * (size - 2.0);
    x_variance_term += size * (size - 1.0) * (2.0 * size + 5.0);
  }
  for (const auto& group_size : tie_group_sizes(y)) {
    double size = static_cast<double>(group_size);
    y_pair_term += size * (size - 1.0);
    y_triplet_term += size * (size - 1.0) * (size - 2.0);
    y_variance_term += size * (size - 1.0) * (2.0 * size + 5.0);
  }
  return (sample_size * (sample_size - 1.0) * (2.0 * sample_size + 5.0) -
          x_variance_term - y_variance_term) /
           18.0 +
         x_pair_term * y_pair_term / (2.0 * sample_size * (sample_size - 1.0)) +
         x_triplet_term * y_triplet_term /
           (9.0 * sample_size * (sample_size - 1.0) * (sample_size - 2.0));
}

void
test_kendall_tie_adjustment()
{
  std::vector<double> x{ 1, 1, 2, 3, 3, 3, 4, 5, 5, 6 };
  std::vector<double> y{ 2, 1, 2, 4, 3, 4, 6, 5, 6, 5 };
  wdm::Indep_test result(x, y, "kendall");
  test::check_near(result.statistic(),
                   kendall_score(x, y) /
                     std::sqrt(kendall_score_variance(x, y)),
                   "Kendall statistic uses the independent tie variance");

  std::vector<double> triple_x{ 1, 1, 1, 2, 3, 4, 5, 6 };
  std::vector<double> triple_y{ 1, 2, 3, 3, 3, 4, 5, 6 };
  wdm::Indep_test triple_result(triple_x, triple_y, "kendall");
  test::check_near(
    triple_result.statistic(),
    kendall_score(triple_x, triple_y) /
      std::sqrt(kendall_score_variance(triple_x, triple_y)),
    "Kendall statistic includes triplet terms from both margins");
}

void
test_tied_triplet_counts()
{
  std::vector<double> values{ 1, 1, 1, 2, 3, 3, 3, 3, 4, 5, 5, 5 };
  std::vector<double> weights{ 1, 2, 3, 9, 1, 2, 3, 4, 8, 2, 3, 5 };
  test::check_near(wdm::utils::count_tied_triplets(values, {}),
                   6.0,
                   "unweighted triplet count spans several tie groups");
  test::check_near(wdm::utils::count_tied_triplets(values, weights),
                   86.0,
                   "weighted triplet count spans several tie groups");
}

void
test_effective_sample_size()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6 };
  std::vector<double> y{ 2, 1, 4, 3, 6, 5 };
  std::vector<double> weights{ 1, 2, 3, 4, 5, 6 };
  wdm::Indep_test weighted(x, y, "pearson", weights);
  test::check_near(
    weighted.n_eff(), 21.0 * 21.0 / 91.0, "Kish effective sample size");
  test::check_near(weighted.statistic(),
                   std::atanh(weighted.estimate()) *
                     std::sqrt(weighted.n_eff() - 3.0),
                   "weighted Pearson uses effective sample size");
  wdm::Indep_test scaled(x, y, "pearson", { 10, 20, 30, 40, 50, 60 });
  test::check_near(scaled.n_eff(),
                   weighted.n_eff(),
                   "effective sample size ignores weight scale");
  wdm::Indep_test uniform(x, y, "pearson", std::vector<double>(x.size(), 3.0));
  test::check_near(uniform.n_eff(),
                   static_cast<double>(x.size()),
                   "uniform weights retain the sample size");

  double missing = std::numeric_limits<double>::quiet_NaN();
  wdm::Indep_test incomplete(
    { 1, 2, 3, 4, 5, 6 }, { 2, missing, 4, 3, 6, 5 }, "pearson", weights);
  test::check_near(incomplete.n_eff(),
                   19.0 * 19.0 / 87.0,
                   "effective sample size follows missing-value removal");

  wdm::Indep_test unweighted_kendall(x, y, "kendall");
  wdm::Indep_test uniform_kendall(
    x, y, "kendall", std::vector<double>(x.size(), 3.0));
  test::check_near(uniform_kendall.statistic(),
                   unweighted_kendall.statistic(),
                   "uniform weights retain the Kendall statistic");
  wdm::Indep_test weighted_kendall(x, y, "kendall", weights);
  wdm::Indep_test scaled_kendall(x, y, "kendall", { 10, 20, 30, 40, 50, 60 });
  test::check_near(scaled_kendall.statistic(),
                   weighted_kendall.statistic(),
                   "Kendall statistic ignores weight scale");
  test::check_near(scaled_kendall.p_value(),
                   weighted_kendall.p_value(),
                   "Kendall p-value ignores weight scale");
  wdm::Indep_test zero_weight_kendall({ 1, 2, 3, 3.5, 4, 5, 6 },
                                      { 2, 1, 4, 100, 3, 6, 5 },
                                      "kendall",
                                      { 1, 2, 3, 0, 4, 5, 6 });
  test::check_near(zero_weight_kendall.statistic(),
                   weighted_kendall.statistic(),
                   "Kendall inference ignores zero-weight rows");
}

void
test_chatterjee_inference_paths()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6 };
  std::vector<double> y{ 2, 5, 1, 6, 3, 4 };
  std::vector<double> weights{ 1, 2, 3, 4, 5, 6 };
  std::tuple<double, double, double, double> continuous =
    wdm::impl::cxi(x, y, weights, true);
  wdm::Indep_test continuous_test(x, y, "chatterjee", weights);
  test::check_near(continuous_test.estimate(),
                   std::get<0>(continuous),
                   "continuous Chatterjee estimate");
  test::check_near(continuous_test.statistic(),
                   (std::get<3>(continuous) - std::get<2>(continuous)) /
                     std::get<1>(continuous),
                   "continuous Chatterjee statistic");

  std::tuple<double, double, double, double> discrete =
    wdm::impl::cxi(x, y, {}, true, "max", {}, false);
  wdm::Indep_test discrete_test(
    x, y, "chatterjee", {}, true, "greater", {}, false);
  test::check_near(discrete_test.statistic(),
                   (std::get<3>(discrete) - std::get<2>(discrete)) /
                     std::get<1>(discrete),
                   "declared-discrete Chatterjee statistic");
  test::check_throws(
    [&]() {
      wdm::Indep_test unavailable(
        x, y, "chatterjee", weights, true, "greater", {}, false);
    },
    "unequally weighted discrete Chatterjee inference is unavailable");

  std::vector<double> tied_y{ 1, 2, 2, 3, 3, 1 };
  wdm::Indep_test tied(x, tied_y, "chatterjee");
  wdm::Indep_test uniformly_weighted(
    x, tied_y, "chatterjee", std::vector<double>(x.size(), 4.0));
  test::check_near(uniformly_weighted.statistic(),
                   tied.statistic(),
                   "uniform weights retain tied-response inference");
  test::check_throws(
    [&]() {
      wdm::Indep_test unavailable(
        x, tied_y, "chatterjee", { 1, 2, 1, 2, 1, 2 });
    },
    "unequal weights reject observed Chatterjee response ties");

  std::vector<double> tied_x{ 1, 1, 2, 2, 3, 3 };
  wdm::Indep_test seeded_first(
    tied_x, y, "chatterjee", {}, true, "greater", { 17, 29 });
  wdm::Indep_test seeded_second(
    tied_x, y, "chatterjee", {}, true, "greater", { 17, 29 });
  test::check_near(seeded_first.statistic(),
                   seeded_second.statistic(),
                   "fixed predictor-tie seeds reproduce Chatterjee inference");
}

void
test_fixed_seed_null_simulation()
{
  const size_t sample_size = 80;
  const size_t replications = 600;
  std::mt19937 generator(20260821);
  std::normal_distribution<double> normal;
  std::vector<std::string> methods{
    "pearson", "spearman", "kendall", "blomqvist", "chatterjee"
  };
  std::vector<double> statistic_sums(methods.size(), 0.0);
  std::vector<double> squared_statistic_sums(methods.size(), 0.0);
  double hoeffding_p_value_sum = 0.0;
  size_t hoeffding_rejections = 0;
  std::vector<double> x(sample_size);
  std::vector<double> y(sample_size);
  for (size_t replication = 0; replication < replications; ++replication) {
    for (size_t i = 0; i < sample_size; ++i) {
      x[i] = normal(generator);
      y[i] = normal(generator);
    }
    for (size_t method_index = 0; method_index < methods.size();
         ++method_index) {
      double statistic =
        wdm::Indep_test(x, y, methods[method_index]).statistic();
      statistic_sums[method_index] += statistic;
      squared_statistic_sums[method_index] += statistic * statistic;
    }
    double hoeffding_p_value = wdm::Indep_test(x, y, "hoeffding").p_value();
    hoeffding_p_value_sum += hoeffding_p_value;
    if (hoeffding_p_value < 0.05)
      ++hoeffding_rejections;
  }

  for (size_t method_index = 0; method_index < methods.size(); ++method_index) {
    double statistic_mean = statistic_sums[method_index] / replications;
    double statistic_variance =
      squared_statistic_sums[method_index] / replications -
      statistic_mean * statistic_mean;
    test::check(std::fabs(statistic_mean) < 0.15,
                methods[method_index] + " fixed-seed null mean");
    test::check(statistic_variance > 0.7 && statistic_variance < 1.3,
                methods[method_index] + " fixed-seed null variance");
  }
  test::check(std::fabs(hoeffding_p_value_sum / replications - 0.5) < 0.2,
              "Hoeffding fixed-seed null p-value mean");
  test::check(static_cast<double>(hoeffding_rejections) / replications < 0.12,
              "Hoeffding fixed-seed null rejection rate");
}

} // namespace

int
main()
{
  test_statistic_transformations();
  test_hoeffding_tail_table();
  test_alternatives_and_p_values();
  test_aliases();
  test_kendall_tie_adjustment();
  test_tied_triplet_counts();
  test_effective_sample_size();
  test_chatterjee_inference_paths();
  test_fixed_seed_null_simulation();
  return test::finish();
}
