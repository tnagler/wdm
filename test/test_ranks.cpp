#include <algorithm>
#include <cmath>
#include <vector>
#include <wdm/ranks.hpp>

#include "test_helpers.hpp"

namespace {

std::vector<double>
sorted(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  return values;
}

std::vector<double>
successive_differences(const std::vector<double>& values)
{
  std::vector<double> differences(values.size());
  for (size_t i = 0; i < values.size(); ++i)
    differences[i] = values[i] - (i > 0 ? values[i - 1] : 0.0);
  return differences;
}

void
test_rank_ties()
{
  std::vector<int> seeds{ 1, 2, 3 };
  std::vector<double> values{ 1, 2, 2, 2, 2, 2, 2, 3, 4, 5 };
  std::vector<double> consecutive{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

  test::check_vector_near(wdm::impl::rank(values, {}, "first"),
                          consecutive,
                          "'first' ranks ties in input order");
  test::check_vector_near(wdm::impl::rank(values, {}, "min"),
                          { 1, 2, 2, 2, 2, 2, 2, 8, 9, 10 },
                          "'min' assigns the minimum rank");
  test::check_vector_near(wdm::impl::rank(values, {}, "average"),
                          { 1, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 8, 9, 10 },
                          "'average' assigns the average rank");

  auto random_ranks = wdm::impl::rank(values, {}, "random", seeds);
  test::check_vector_near(
    sorted(random_ranks), consecutive, "'random' permutes consecutive ranks");
  test::check_vector_near(random_ranks,
                          wdm::impl::rank(values, {}, "random", seeds),
                          "'random' is reproducible with fixed seeds");

  std::vector<double> tied(5, 2.0);
  std::vector<double> weights{ 1, 2, 3, 4, 5 };
  auto weighted = sorted(wdm::impl::rank(tied, weights, "random", seeds));
  test::check_vector_near(
    sorted(successive_differences(weighted)),
    { 1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0 },
    "weighted random ranks use cumulative normalized weights");
  test::check_near(weighted.back(), 5.0, "largest weighted rank");

  std::vector<double> pairs;
  for (size_t i = 0; i < 20; ++i)
    pairs.insert(pairs.end(), 2, static_cast<double>(i));
  auto shuffled = wdm::impl::rank(pairs, {}, "random", seeds);
  bool ascending = false;
  bool descending = false;
  for (size_t i = 0; i < shuffled.size(); i += 2) {
    ascending = ascending || shuffled[i] < shuffled[i + 1];
    descending = descending || shuffled[i] > shuffled[i + 1];
  }
  test::check(ascending && descending, "tie groups are shuffled independently");
  test::check_vector_near(sorted(successive_differences(sorted(shuffled))),
                          std::vector<double>(pairs.size(), 1.0),
                          "random tie breaking removes ties");

  std::vector<double> with_nan{ 2, NAN, 2, 2 };
  auto nan_ranks = wdm::impl::rank(with_nan, {}, "random", seeds);
  test::check(std::isnan(nan_ranks[1]), "random ranks preserve NaNs");
  nan_ranks.erase(nan_ranks.begin() + 1);
  test::check_vector_near(
    sorted(nan_ranks), { 1, 2, 3 }, "random ranks ignore NaNs");
}

void
test_rank0_ties()
{
  std::vector<double> values{ 1, 3, 2, 5, 3, 2, 20, 15 };
  test::check_vector_near(
    wdm::impl::rank0(values), { 0, 3, 1, 5, 3, 1, 7, 6 }, "rank0 'min'");
  test::check_vector_near(wdm::impl::rank0(values, {}, "max"),
                          { 1, 5, 3, 6, 5, 3, 8, 7 },
                          "rank0 'max'");
  test::check_throws([&]() { wdm::impl::rank0(values, {}, "unknown"); },
                     "rank0 rejects an unknown ties method");

  values = { 1, 2, 2, 2, 2, 2, 2, 3, 4, 5 };
  std::vector<double> weights{ 1, 1, 2, 2, 1, 3, 1, 1, 1, 1 };
  test::check_vector_near(wdm::impl::rank0(values, weights),
                          { 0, 1, 1, 1, 1, 1, 1, 11, 12, 13 },
                          "weighted rank0 accumulates weights");
  test::check_vector_near(wdm::impl::rank0(values, {}, "average"),
                          { 0, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 7, 8, 9 },
                          "unweighted rank0 average ties");
  test::check_vector_near(wdm::impl::rank0(values, weights, "average"),
                          { 0, 5, 5, 5, 5, 5, 5, 11, 12, 13 },
                          "weighted rank0 average ties");
}

}

int
main()
{
  test_rank_ties();
  test_rank0_ties();
  return test::finish();
}
