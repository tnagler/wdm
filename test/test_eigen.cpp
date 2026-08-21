#include <cmath>
#include <iostream>
#include <wdm/eigen.hpp>

int
main()
{
  Eigen::MatrixXd observations(40, 2);
  size_t row = 0;
  for (int value = -20; value <= 20; ++value) {
    if (value == 0)
      continue;
    observations(row, 0) = value;
    observations(row, 1) = value * value;
    ++row;
  }

  std::vector<int> seeds{ 17, 29, 43 };
  Eigen::MatrixXd chatterjee =
    wdm::wdm(observations, "cxi", Eigen::VectorXd(), true, seeds);
  double forward =
    wdm::wdm(observations.col(0), observations.col(1), "cxi", {}, true, seeds);
  double reverse =
    wdm::wdm(observations.col(1), observations.col(0), "cxi", {}, true, seeds);
  if (std::fabs(chatterjee(0, 1) - forward) > 1e-12 ||
      std::fabs(chatterjee(1, 0) - reverse) > 1e-12 ||
      chatterjee(0, 1) <= chatterjee(1, 0) + 0.5) {
    std::cout << "FAILED: Chatterjee matrix preserves both directions"
              << std::endl;
    return 1;
  }

  Eigen::MatrixXd pearson = wdm::wdm(observations, "pearson");
  if (pearson(0, 1) != pearson(1, 0)) {
    std::cout << "FAILED: non-Chatterjee matrix remains symmetric" << std::endl;
    return 1;
  }

  std::cout << "all Eigen checks passed" << std::endl;
  return 0;
}
