#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace test {

static int failures = 0;

inline void
check(bool condition, const std::string& description)
{
  if (!condition) {
    std::cout << "FAILED: " << description << std::endl;
    ++failures;
  }
}

inline void
check_near(double actual,
           double expected,
           const std::string& description,
           double tolerance = 1e-10)
{
  if (!(std::fabs(actual - expected) <= tolerance)) {
    std::cout << "FAILED: " << description << " (got " << std::setprecision(12)
              << actual << ", expected " << expected << ")" << std::endl;
    ++failures;
  }
}

inline void
check_vector_near(const std::vector<double>& actual,
                  const std::vector<double>& expected,
                  const std::string& description,
                  double tolerance = 1e-10)
{
  if (actual.size() != expected.size()) {
    check(false, description + " (different sizes)");
    return;
  }
  for (size_t i = 0; i < actual.size(); ++i) {
    if (!(std::fabs(actual[i] - expected[i]) <= tolerance)) {
      check_near(actual[i],
                 expected[i],
                 description + " at index " + std::to_string(i),
                 tolerance);
      return;
    }
  }
}

template<typename Function>
void
check_throws(Function function, const std::string& description)
{
  bool threw = false;
  try {
    function();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, description);
}

inline int
finish()
{
  if (failures > 0) {
    std::cout << failures << " check(s) failed" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}

}
