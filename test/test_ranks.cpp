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

//! 13 tie groups of about 150 values each: large enough that an unstable sort
//! reorders the values within a group.
std::vector<double>
many_ties()
{
  std::vector<double> x;
  for (int i = 0; i < 2000; i++)
    x.push_back(static_cast<double>((i * 7919) % 13));
  return x;
}

void
test_order_within_ties()
{
  auto x = many_ties();
  auto order = wdm::utils::get_order(x);
  bool stable = true;
  for (size_t i = 1; i < order.size(); ++i)
    if (x[order[i]] == x[order[i - 1]])
      stable = stable && order[i] > order[i - 1];
  test::check(stable, "get_order keeps equal values in order of appearance");

  auto first = wdm::impl::rank(x, {}, "first");
  bool in_input_order = true;
  for (size_t i = 1; i < x.size(); ++i)
    for (size_t j = 0; j < i && in_input_order; ++j)
      if (x[j] == x[i])
        in_input_order = first[j] < first[i];
  test::check(in_input_order, "'first' ranks many ties in input order");
}

void
test_tie_order()
{
  auto x = many_ties();
  std::vector<int> seeds{ 5 };
  auto order = wdm::utils::get_order(x);

  std::vector<size_t> sizes;
  for (size_t i = 0; i < order.size();) {
    size_t m = 1;
    while (i + m < order.size() && x[order[i + m]] == x[order[i]])
      ++m;
    sizes.push_back(m);
    i += m;
  }

  auto ord = wdm::impl::tie_order(sizes, seeds);
  std::vector<double> ranks(x.size());
  for (size_t g = 0, o = 0; g < sizes.size(); o += sizes[g++])
    for (size_t k = 0; k < sizes[g]; ++k)
      ranks[order[o + ord[o + k]]] = static_cast<double>(o + k + 1);
  test::check_vector_near(ranks,
                          wdm::impl::rank(x, {}, "random", seeds),
                          "tie_order reproduces the random ranks");

  test::check_vector_near(
    [&] {
      std::vector<double> v;
      for (size_t i : wdm::impl::tie_order({ 1, 1, 1 }, seeds))
        v.push_back(static_cast<double>(i));
      return v;
    }(),
    { 0, 0, 0 },
    "values without ties keep their place");
  test::check_throws([&]() { wdm::impl::tie_order({ 2, 0 }, seeds); },
                     "tie_order refuses an empty group");
}

#ifdef USE_BOOST
void
test_random_ranks_are_portable()
{
  // The Boost generator and distributions are specified exactly, so with a
  // stable order within ties the random ranks are the same on every platform.
  std::vector<double> x;
  for (int i = 0; i < 40; i++)
    x.push_back(static_cast<double>((i * 7) % 4));
  test::check_vector_near(wdm::impl::rank(x, {}, "random", { 5 }),
                          { 6,  33, 23, 16, 1,  36, 24, 14, 3,  32,
                            21, 18, 4,  37, 28, 12, 2,  38, 30, 19,
                            7,  31, 27, 13, 8,  40, 22, 20, 9,  39,
                            25, 11, 10, 35, 29, 17, 5,  34, 26, 15 },
                          "random ranks match the pinned values");
}
#endif

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
  test_order_within_ties();
  test_tie_order();
#ifdef USE_BOOST
  test_random_ranks_are_portable();
#endif
  test_rank0_ties();
  return test::finish();
}
