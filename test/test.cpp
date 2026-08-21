/*
int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/

// #include <wdm/include/wdm.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
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

  wdm::Indep_test test(v, v_sq, "cxi");
  check(std::isfinite(test.p_value()), "xi p-value is finite");
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

  test_rank0_ties_methods();
  test_cxi();

  if (failures > 0) {
    std::cout << failures << " check(s) failed" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}
