#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>
#include <wdm.hpp>

#include "test_helpers.hpp"

namespace {

std::vector<double>
case_weights(size_t sample_size, std::vector<double> weights)
{
  if (weights.empty())
    weights.assign(sample_size, 1.0);
  return weights;
}

double
reference_pearson(const std::vector<double>& x,
                  const std::vector<double>& y,
                  std::vector<double> weights)
{
  weights = case_weights(x.size(), weights);
  double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  double x_mean = 0.0;
  double y_mean = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    x_mean += weights[i] * x[i] / weight_sum;
    y_mean += weights[i] * y[i] / weight_sum;
  }

  double covariance = 0.0;
  double x_variance = 0.0;
  double y_variance = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    covariance += weights[i] * (x[i] - x_mean) * (y[i] - y_mean);
    x_variance += weights[i] * (x[i] - x_mean) * (x[i] - x_mean);
    y_variance += weights[i] * (y[i] - y_mean) * (y[i] - y_mean);
  }
  return covariance / std::sqrt(x_variance * y_variance);
}

std::vector<double>
reference_average_ranks(const std::vector<double>& values,
                        std::vector<double> weights)
{
  weights = case_weights(values.size(), weights);
  std::vector<double> ranks(values.size(), 0.0);
  for (size_t i = 0; i < values.size(); ++i) {
    double tied_weight = 0.0;
    double tied_pair_weight = 0.0;
    for (size_t j = 0; j < values.size(); ++j) {
      if (values[j] < values[i])
        ranks[i] += weights[j];
      if (values[j] == values[i]) {
        tied_pair_weight += tied_weight * weights[j];
        tied_weight += weights[j];
      }
    }
    ranks[i] += tied_pair_weight / tied_weight;
  }
  return ranks;
}

double
reference_spearman(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& weights)
{
  return reference_pearson(reference_average_ranks(x, weights),
                           reference_average_ranks(y, weights),
                           weights);
}

double
reference_kendall(const std::vector<double>& x,
                  const std::vector<double>& y,
                  std::vector<double> weights)
{
  weights = case_weights(x.size(), weights);
  double concordance_difference = 0.0;
  double untied_x_pair_weight = 0.0;
  double untied_y_pair_weight = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = i + 1; j < x.size(); ++j) {
      double pair_weight = weights[i] * weights[j];
      if (x[i] != x[j])
        untied_x_pair_weight += pair_weight;
      if (y[i] != y[j])
        untied_y_pair_weight += pair_weight;
      if (x[i] != x[j] && y[i] != y[j])
        concordance_difference +=
          pair_weight * ((x[i] < x[j]) == (y[i] < y[j]) ? 1.0 : -1.0);
    }
  }
  return concordance_difference /
         std::sqrt(untied_x_pair_weight * untied_y_pair_weight);
}

double
reference_median(const std::vector<double>& values, std::vector<double> weights)
{
  weights = case_weights(values.size(), weights);
  std::vector<size_t> order(values.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
    return values[i] < values[j];
  });

  double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
  double pair_weight = 0.0;
  double preceding_weight = 0.0;
  for (size_t i = 0; i < weights.size(); ++i) {
    pair_weight += preceding_weight * weights[i];
    preceding_weight += weights[i];
  }
  double median_rank = pair_weight / total_weight;
  preceding_weight = 0.0;
  for (size_t position = 0; position < order.size(); ++position) {
    if (preceding_weight == median_rank)
      return values[order[position]];
    if (preceding_weight > median_rank)
      return 0.5 * (values[order[position - 1]] + values[order[position]]);
    preceding_weight += weights[order[position]];
  }
  return values[order.back()];
}

double
reference_blomqvist(const std::vector<double>& x,
                    const std::vector<double>& y,
                    std::vector<double> weights)
{
  double x_median = reference_median(x, weights);
  double y_median = reference_median(y, weights);
  weights = case_weights(x.size(), weights);
  double matching_quadrant_weight = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    if ((x[i] <= x_median && y[i] <= y_median) ||
        (x[i] > x_median && y[i] > y_median))
      matching_quadrant_weight += weights[i];
  }
  return 2.0 * matching_quadrant_weight /
           std::accumulate(weights.begin(), weights.end(), 0.0) -
         1.0;
}

double
ordered_weight_product_sum(const std::vector<double>& weights,
                           size_t factor_count,
                           size_t first_index = 0)
{
  if (factor_count == 0)
    return 1.0;
  double product_sum = 0.0;
  for (size_t i = first_index; i + factor_count <= weights.size(); ++i)
    product_sum +=
      weights[i] * ordered_weight_product_sum(weights, factor_count - 1, i + 1);
  return static_cast<double>(factor_count) * product_sum;
}

std::vector<double>
reference_rank_below(const std::vector<double>& values,
                     const std::vector<double>& weights,
                     unsigned int weight_power)
{
  std::vector<double> ranks(values.size(), 0.0);
  for (size_t i = 0; i < values.size(); ++i) {
    for (size_t j = 0; j < values.size(); ++j) {
      if (values[j] < values[i])
        ranks[i] += std::pow(weights[j], static_cast<int>(weight_power));
    }
  }
  return ranks;
}

std::vector<double>
reference_bivariate_rank(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& weights,
                         unsigned int weight_power)
{
  std::vector<double> ranks(x.size(), 0.0);
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < x.size(); ++j) {
      if (x[j] < x[i] && y[j] < y[i])
        ranks[i] += std::pow(weights[j], static_cast<int>(weight_power));
    }
  }
  return ranks;
}

double
reference_hoeffding(const std::vector<double>& x,
                    const std::vector<double>& y,
                    std::vector<double> weights)
{
  weights = case_weights(x.size(), weights);
  std::vector<double> x_rank = reference_rank_below(x, weights, 1);
  std::vector<double> y_rank = reference_rank_below(y, weights, 1);
  std::vector<double> squared_x_rank = reference_rank_below(x, weights, 2);
  std::vector<double> squared_y_rank = reference_rank_below(y, weights, 2);
  std::vector<double> joint_rank = reference_bivariate_rank(x, y, weights, 1);
  std::vector<double> squared_joint_rank =
    reference_bivariate_rank(x, y, weights, 2);
  std::vector<double> cubed_joint_rank =
    reference_bivariate_rank(x, y, weights, 3);
  std::vector<double> fourth_joint_rank =
    reference_bivariate_rank(x, y, weights, 4);

  double first_component = 0.0;
  double second_component = 0.0;
  double third_component = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    first_component +=
      weights[i] * (joint_rank[i] * joint_rank[i] - squared_joint_rank[i]);
    second_component +=
      weights[i] *
      ((x_rank[i] * y_rank[i] - squared_joint_rank[i]) * joint_rank[i] -
       squared_joint_rank[i] * (x_rank[i] + y_rank[i]) +
       2.0 * cubed_joint_rank[i]);
    third_component +=
      weights[i] * ((x_rank[i] * x_rank[i] - squared_x_rank[i]) *
                      (y_rank[i] * y_rank[i] - squared_y_rank[i]) -
                    4.0 * ((x_rank[i] * y_rank[i] - squared_joint_rank[i]) *
                             squared_joint_rank[i] -
                           cubed_joint_rank[i] * (x_rank[i] + y_rank[i]) +
                           2.0 * fourth_joint_rank[i]) -
                    2.0 * (squared_joint_rank[i] * squared_joint_rank[i] -
                           fourth_joint_rank[i]));
  }
  return 30.0 *
         (first_component / ordered_weight_product_sum(weights, 3) -
          2.0 * second_component / ordered_weight_product_sum(weights, 4) +
          third_component / ordered_weight_product_sum(weights, 5));
}

double
reference_chatterjee(std::vector<double> x,
                     std::vector<double> y,
                     std::vector<double> weights)
{
  weights = case_weights(x.size(), weights);
  std::vector<size_t> order(x.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
    return x[i] < x[j];
  });
  double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::vector<double> distribution_rank(x.size(), 0.0);
  std::vector<double> survival_rank(x.size(), 0.0);
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < x.size(); ++j) {
      if (y[j] <= y[i])
        distribution_rank[i] += weights[j] / weight_sum;
      if (y[j] >= y[i])
        survival_rank[i] += weights[j] / weight_sum;
    }
  }

  double edge_difference = 0.0;
  for (size_t position = 0; position + 1 < order.size(); ++position) {
    edge_difference += weights[order[position]] / weight_sum *
                       std::fabs(distribution_rank[order[position + 1]] -
                                 distribution_rank[order[position]]);
  }
  double denominator = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    denominator += 2.0 * weights[i] / weight_sum * survival_rank[i] *
                   (1.0 - survival_rank[i]);
  }
  return 1.0 - edge_difference / denominator;
}

void
test_exact_references()
{
  std::vector<double> x{ 4, 1, 7, 2, 6, 3, 5 };
  std::vector<double> y{ 5, 2, 1, 6, 3, 7, 4 };
  std::vector<double> weights{ 1, 2, 1, 3, 2, 4, 1 };

  test::check_near(wdm::wdm(x, y, "pearson"),
                   reference_pearson(x, y, {}),
                   "unweighted Pearson reference");
  test::check_near(wdm::wdm(x, y, "pearson", weights),
                   reference_pearson(x, y, weights),
                   "weighted Pearson reference");
  test::check_near(wdm::wdm(x, y, "spearman"),
                   reference_spearman(x, y, {}),
                   "unweighted Spearman reference");
  test::check_near(wdm::wdm(x, y, "spearman", weights),
                   reference_spearman(x, y, weights),
                   "weighted Spearman reference");
  test::check_near(wdm::wdm(x, y, "kendall"),
                   reference_kendall(x, y, {}),
                   "unweighted Kendall reference");
  test::check_near(wdm::wdm(x, y, "kendall", weights),
                   reference_kendall(x, y, weights),
                   "weighted Kendall reference");
  test::check_near(wdm::wdm(x, y, "blomqvist"),
                   reference_blomqvist(x, y, {}),
                   "unweighted Blomqvist reference");
  test::check_near(wdm::wdm(x, y, "blomqvist", weights),
                   reference_blomqvist(x, y, weights),
                   "weighted Blomqvist reference");
  test::check_near(wdm::wdm(x, y, "hoeffding"),
                   reference_hoeffding(x, y, {}),
                   "unweighted Hoeffding reference");
  test::check_near(wdm::wdm(x, y, "hoeffding", weights),
                   reference_hoeffding(x, y, weights),
                   "weighted Hoeffding reference");
  test::check_near(wdm::wdm(x, y, "chatterjee"),
                   reference_chatterjee(x, y, {}),
                   "unweighted Chatterjee reference");
  test::check_near(wdm::wdm(x, y, "chatterjee", weights),
                   reference_chatterjee(x, y, weights),
                   "weighted Chatterjee reference");
}

void
test_aliases()
{
  std::vector<double> x{ 4, 1, 7, 2, 6, 3, 5 };
  std::vector<double> y{ 5, 2, 1, 6, 3, 7, 4 };
  std::vector<std::vector<std::string>> aliases{
    { "pearson", "prho", "cor" },   { "spearman", "srho", "rho" },
    { "kendall", "ktau", "tau" },   { "blomqvist", "bbeta", "beta" },
    { "hoeffding", "hoeffd", "d" }, { "chatterjee", "cxi", "xi" }
  };
  for (const auto& method_aliases : aliases) {
    double expected = wdm::wdm(x, y, method_aliases.front());
    for (const auto& method : method_aliases)
      test::check_near(wdm::wdm(x, y, method), expected, method + " alias");
  }
}

void
test_invariances()
{
  std::vector<double> x{ 4, 1, 7, 2, 6, 3, 5 };
  std::vector<double> y{ 5, 2, 1, 6, 3, 7, 4 };
  std::vector<double> weights{ 1, 2, 1, 3, 2, 4, 1 };
  std::vector<double> scaled_weights{ 10, 20, 10, 30, 20, 40, 10 };
  std::vector<double> permuted_x{ 3, 6, 1, 5, 4, 7, 2 };
  std::vector<double> permuted_y{ 7, 3, 2, 4, 5, 1, 6 };
  std::vector<double> permuted_weights{ 4, 2, 2, 1, 1, 1, 3 };
  std::vector<std::string> methods{ "pearson",   "spearman",  "kendall",
                                    "blomqvist", "hoeffding", "chatterjee" };

  for (const auto& method : methods) {
    double unweighted = wdm::wdm(x, y, method);
    double weighted = wdm::wdm(x, y, method, weights);
    test::check_near(wdm::wdm(x, y, method, std::vector<double>(x.size(), 1.0)),
                     unweighted,
                     method + " uniform weights");
    test::check_near(wdm::wdm(x, y, method, scaled_weights),
                     weighted,
                     method + " weight scaling");
    test::check_near(wdm::wdm(permuted_x, permuted_y, method, permuted_weights),
                     weighted,
                     method + " row permutation");
    if (method != "chatterjee")
      test::check_near(
        wdm::wdm(y, x, method, weights), weighted, method + " symmetry");
  }
}

void
test_zero_weight_equivalence()
{
  std::vector<double> x{ 4, 1, 7, 2, 6, 3, 5 };
  std::vector<double> y{ 5, 2, 1, 6, 3, 7, 4 };
  std::vector<double> weights{ 1, 2, 1, 3, 2, 4, 1 };
  std::vector<double> augmented_x{ 4, 1, 7, 2, 6, 3, 5, 4.5 };
  std::vector<double> augmented_y{ 5, 2, 1, 6, 3, 7, 4, 100 };
  std::vector<double> augmented_weights{ 1, 2, 1, 3, 2, 4, 1, 0 };
  std::vector<std::string> methods{ "pearson",   "spearman",  "kendall",
                                    "blomqvist", "hoeffding", "chatterjee" };
  for (const auto& method : methods) {
    test::check_near(
      wdm::wdm(augmented_x, augmented_y, method, augmented_weights),
      wdm::wdm(x, y, method, weights),
      method + " ignores zero-weight rows");
  }
}

void
test_ties_and_endpoints()
{
  std::vector<double> tied_x{ 1, 1, 2, 3, 3, 4, 5, 5 };
  std::vector<double> tied_y{ 3, 2, 2, 4, 1, 4, 5, 5 };
  std::vector<double> weights{ 1, 2, 3, 1, 2, 4, 1, 2 };
  test::check_near(wdm::wdm(tied_x, tied_y, "spearman"),
                   reference_spearman(tied_x, tied_y, {}),
                   "Spearman average response ties");
  test::check_near(wdm::wdm(tied_x, tied_y, "spearman", weights),
                   reference_spearman(tied_x, tied_y, weights),
                   "weighted Spearman average ties");
  test::check_near(wdm::wdm(tied_x, tied_y, "kendall"),
                   reference_kendall(tied_x, tied_y, {}),
                   "Kendall tau-b ties");
  test::check_near(wdm::wdm(tied_x, tied_y, "kendall", weights),
                   reference_kendall(tied_x, tied_y, weights),
                   "weighted Kendall tau-b ties");
  test::check_near(wdm::wdm(tied_x, tied_y, "blomqvist"),
                   reference_blomqvist(tied_x, tied_y, {}),
                   "Blomqvist median ties");

  std::vector<double> increasing{ 1, 2, 3, 4, 5, 6, 7, 8 };
  std::vector<double> decreasing(increasing.rbegin(), increasing.rend());
  for (const auto& method : std::vector<std::string>{
         "pearson", "spearman", "kendall", "blomqvist" }) {
    test::check_near(wdm::wdm(increasing, increasing, method),
                     1.0,
                     method + " upper endpoint");
    test::check_near(wdm::wdm(increasing, decreasing, method),
                     -1.0,
                     method + " lower endpoint");
  }
  test::check_near(wdm::wdm(increasing, increasing, "hoeffding"),
                   1.0,
                   "Hoeffding upper endpoint");
  test::check_near(wdm::wdm(increasing, decreasing, "hoeffding"),
                   1.0,
                   "Hoeffding detects decreasing functional dependence");
  test::check_near(wdm::wdm(increasing, increasing, "chatterjee"),
                   2.0 / 3.0,
                   "finite-sample Chatterjee increasing endpoint");
  test::check_near(wdm::wdm(increasing, decreasing, "chatterjee"),
                   2.0 / 3.0,
                   "finite-sample Chatterjee decreasing endpoint");

  std::vector<double> nonlinear_x;
  std::vector<double> nonlinear_y;
  for (int value = -20; value <= 20; ++value) {
    if (value != 0) {
      nonlinear_x.push_back(value);
      nonlinear_y.push_back(value * value);
    }
  }
  test::check(wdm::wdm(nonlinear_x, nonlinear_y, "chatterjee") >
                wdm::wdm(nonlinear_y, nonlinear_x, "chatterjee") + 0.5,
              "Chatterjee is asymmetric");
}

} // namespace

int
main()
{
  test_exact_references();
  test_aliases();
  test_invariances();
  test_zero_weight_equivalence();
  test_ties_and_endpoints();
  return test::finish();
}
