#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <wdm.hpp>

#include "test_helpers.hpp"

namespace {

const double missing = std::numeric_limits<double>::quiet_NaN();

void
test_size_mismatches()
{
  test::check_throws([]() { wdm::wdm({ 1, 2, 3 }, { 1, 2 }, "pearson"); },
                     "estimate rejects different x and y sizes");
  test::check_throws(
    []() { wdm::wdm({ 1, 2, 3 }, { 1, 2, 3 }, "pearson", { 1, 2 }); },
    "estimate rejects a different weight size");
  test::check_throws(
    []() { wdm::Indep_test({ 1, 2, 3 }, { 1, 2 }, "pearson"); },
    "inference rejects different x and y sizes");
  test::check_throws(
    []() { wdm::Indep_test({ 1, 2, 3 }, { 1, 2, 3 }, "pearson", { 1, 2 }); },
    "inference rejects a different weight size");
}

void
test_sample_sizes()
{
  std::vector<std::string> ordinary_methods{
    "pearson", "spearman", "kendall", "blomqvist", "chatterjee"
  };
  for (const auto& method : ordinary_methods) {
    test::check(std::isnan(wdm::wdm({}, {}, method)),
                method + " returns NaN for empty input");
    test::check(std::isnan(wdm::wdm({ 1 }, { 2 }, method)),
                method + " returns NaN for one observation");
    test::check_throws([&]() { wdm::wdm({}, {}, method, {}, false); },
                       method +
                         " rejects empty input when removal is disabled");
    test::check_throws([&]() { wdm::wdm({ 1 }, { 2 }, method, {}, false); },
                       method +
                         " rejects one observation when removal is disabled");
  }

  for (const auto& method :
       std::vector<std::string>{ "hoeffding", "hoeffd", "d" }) {
    test::check(std::isnan(wdm::wdm({ 1, 2, 3, 4 }, { 4, 1, 3, 2 }, method)),
                method + " returns NaN below five observations");
    test::check_throws(
      [&]() { wdm::wdm({ 1, 2, 3, 4 }, { 4, 1, 3, 2 }, method, {}, false); },
      method + " rejects fewer than five observations");
  }
}

void
test_missing_values()
{
  std::vector<double> complete_x{ 1, 3, 4, 5, 6 };
  std::vector<double> complete_y{ 6, 2, 3, 5, 1 };
  std::vector<double> complete_weights{ 1, 3, 4, 5, 6 };
  std::vector<double> incomplete_x{ 1, missing, 3, 4, 5, 6 };
  std::vector<double> incomplete_y{ 6, 4, 2, 3, 5, 1 };
  std::vector<double> incomplete_weights{ 1, 2, 3, 4, 5, 6 };

  test::check_near(
    wdm::wdm(incomplete_x, incomplete_y, "pearson", incomplete_weights),
    wdm::wdm(complete_x, complete_y, "pearson", complete_weights),
    "missing x row is removed");
  incomplete_x[1] = 2;
  incomplete_y[1] = missing;
  test::check_near(
    wdm::wdm(incomplete_x, incomplete_y, "pearson", incomplete_weights),
    wdm::wdm(complete_x, complete_y, "pearson", complete_weights),
    "missing y row is removed");
  incomplete_y[1] = 4;
  incomplete_weights[1] = missing;
  test::check_near(
    wdm::wdm(incomplete_x, incomplete_y, "pearson", incomplete_weights),
    wdm::wdm(complete_x, complete_y, "pearson", complete_weights),
    "missing weight row is removed");

  test::check_throws(
    [&]() { wdm::wdm({ 1, 2, 3 }, { 1, missing, 3 }, "pearson", {}, false); },
    "missing response is rejected when removal is disabled");
  test::check_throws(
    [&]() {
      wdm::wdm({ 1, 2, 3 }, { 3, 1, 2 }, "pearson", { 1, missing, 1 }, false);
    },
    "missing weight is rejected when removal is disabled");
}

void
test_all_observations_removed()
{
  std::vector<double> all_missing{
    missing, missing, missing, missing, missing
  };
  test::check(std::isnan(wdm::wdm(all_missing, all_missing, "pearson")),
              "all-missing estimate is NaN");
  test::check(std::isnan(wdm::wdm(all_missing, all_missing, "hoeffding")),
              "all-missing Hoeffding estimate is NaN");

  wdm::Indep_test result(all_missing, all_missing, "pearson");
  test::check_near(result.n_eff(), 0.0, "all-missing effective sample size");
  test::check(std::isnan(result.estimate()),
              "all-missing inference estimate is NaN");
  test::check(std::isnan(result.statistic()),
              "all-missing inference statistic is NaN");
  test::check(std::isnan(result.p_value()),
              "all-missing inference p-value is NaN");
}

void
test_invalid_weights()
{
  std::vector<double> x{ 1, 2, 3, 4, 5, 6 };
  std::vector<double> y{ 6, 2, 4, 1, 5, 3 };
  std::vector<std::string> methods{ "pearson",   "spearman",  "kendall",
                                    "blomqvist", "hoeffding", "chatterjee" };
  for (const auto& method : methods) {
    test::check_throws([&]() { wdm::wdm(x, y, method, { 1, 1, -1, 1, 1, 1 }); },
                       method + " rejects negative weights");
    test::check_throws(
      [&]() {
        wdm::wdm(x,
                 y,
                 method,
                 { 1, 1, std::numeric_limits<double>::infinity(), 1, 1, 1 });
      },
      method + " rejects infinite weights");
    test::check_throws(
      [&]() { wdm::wdm(x, y, method, std::vector<double>(x.size(), 0.0)); },
      method + " rejects zero total weight");
    test::check_throws(
      [&]() { wdm::Indep_test(x, y, method, { 1, 1, -1, 1, 1, 1 }); },
      method + " inference rejects negative weights");
  }
}

void
test_unknown_methods()
{
  test::check_throws([]() { wdm::wdm({ 1, 2, 3 }, { 3, 1, 2 }, "unknown"); },
                     "estimate rejects an unknown method");
  test::check_throws([]() { wdm::wdm({}, {}, "unknown"); },
                     "unknown method is rejected before sample-size handling");
  test::check_throws(
    []() { wdm::Indep_test({ 1, 2, 3 }, { 3, 1, 2 }, "unknown"); },
    "inference rejects an unknown method");
}

void
test_constant_inputs()
{
  std::vector<double> constant(6, 1.0);
  std::vector<double> varying{ 1, 2, 3, 4, 5, 6 };
  for (const auto& method :
       std::vector<std::string>{ "pearson", "spearman", "kendall" }) {
    test::check(std::isnan(wdm::wdm(constant, varying, method)),
                method + " is undefined for a constant margin");
  }
  test::check_near(wdm::wdm(constant, varying, "blomqvist"),
                   0.0,
                   "Blomqvist constant-margin behavior");
  test::check(std::isfinite(wdm::wdm(constant, varying, "hoeffding")),
              "Hoeffding constant-margin behavior is finite");
  test::check_throws([&]() { wdm::wdm(varying, constant, "chatterjee"); },
                     "Chatterjee rejects a constant response");
  test::check(std::isfinite(
                wdm::wdm(constant, varying, "chatterjee", {}, true, { 1, 2 })),
              "Chatterjee permits a constant predictor");

  wdm::Indep_test result(constant, varying, "pearson");
  test::check(std::isnan(result.estimate()),
              "constant Pearson inference estimate is NaN");
  test::check(std::isnan(result.statistic()),
              "constant Pearson inference statistic is NaN");
  test::check(std::isnan(result.p_value()),
              "constant Pearson inference p-value is NaN");
}

} // namespace

int
main()
{
  test_size_mismatches();
  test_sample_sizes();
  test_missing_values();
  test_all_observations_removed();
  test_invalid_weights();
  test_unknown_methods();
  test_constant_inputs();
  return test::finish();
}
