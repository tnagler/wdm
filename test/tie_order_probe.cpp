// Prints fingerprints of the tie-order-dependent results, so CI logs can be
// compared across platforms.
#include <iostream>
#include <vector>
#include <wdm/ranks.hpp>

namespace {

unsigned long long
fingerprint(const std::vector<double>& values)
{
  unsigned long long h = 1469598103934665603ULL;
  for (double v : values)
    h = (h ^ static_cast<unsigned long long>(v)) * 1099511628211ULL;
  return h;
}

}

int
main()
{
  std::vector<double> x;
  for (int i = 0; i < 2000; i++)
    x.push_back(static_cast<double>((i * 7919) % 13));
  std::vector<double> order;
  for (size_t i : wdm::utils::get_order(x))
    order.push_back(static_cast<double>(i));
  std::cout << "get_order  " << fingerprint(order) << "\n";
  std::cout << "first      " << fingerprint(wdm::impl::rank(x, {}, "first"))
            << "\n";
  std::cout << "random {5} "
            << fingerprint(wdm::impl::rank(x, {}, "random", { 5 })) << "\n";
  return 0;
}
