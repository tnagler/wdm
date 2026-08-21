#include <wdm.hpp>

#include <cmath>
#include <vector>

int
main()
{
  std::vector<double> values{ 1, 2, 3, 4 };
  return std::fabs(wdm::wdm(values, values, "pearson") - 1.0) > 1e-12;
}
