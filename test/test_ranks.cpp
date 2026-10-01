#include <algorithm>
#include <cmath>
#include <limits>
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
test_tie_keys()
{
  auto x = many_ties();
  std::vector<int> seeds{ 5 };
  auto keys = wdm::impl::tie_keys(x.size(), seeds);
  auto sorted_keys = keys;
  std::sort(sorted_keys.begin(), sorted_keys.end());
  bool permutation = true;
  for (size_t i = 0; i < sorted_keys.size(); ++i)
    permutation = permutation && (sorted_keys[i] == i);
  test::check(permutation, "tie_keys is a permutation");

  auto ranks = wdm::impl::rank(x, {}, "random", seeds);
  bool by_key = true;
  for (size_t i = 0; i < x.size(); ++i)
    for (size_t j = 0; j < x.size(); ++j)
      if ((x[i] == x[j]) && (keys[i] < keys[j]))
        by_key = by_key && (ranks[i] < ranks[j]);
  test::check(by_key, "random ranks order ties by the keys");
}

void
test_random_ties_are_local()
{
  // A value that joins a tie group moves the ranks of that group only: every
  // other group keeps its order, however the groups are laid out.
  auto x = many_ties();
  std::vector<int> seeds{ 5 };
  auto before = wdm::impl::rank(x, {}, "random", seeds);
  size_t moved = 0;
  while (x[moved] == x[moved + 1])
    ++moved;
  const double from = x[moved], to = x[moved + 1];
  x[moved] = to;
  auto after = wdm::impl::rank(x, {}, "random", seeds);
  bool others_kept = true;
  for (size_t i = 0; i < x.size(); ++i)
    for (size_t j = 0; j < x.size(); ++j)
      if ((i != moved) && (j != moved) && (x[i] == x[j]) && (x[i] != to) &&
          (x[i] != from))
        others_kept =
          others_kept && ((before[i] < before[j]) == (after[i] < after[j]));
  test::check(others_kept, "a changed tie group leaves the others' order");
}

void
test_random_ranks_are_portable()
{
  // The generator, its seeding and its distributions are specified exactly,
  // and the order within ties is total, so the random ranks are the same on
  // every platform.
  std::vector<double> x;
  for (int i = 0; i < 40; i++)
    x.push_back(static_cast<double>((i * 7) % 4));
  test::check_vector_near(wdm::impl::rank(x, {}, "random", { 5 }),
                          { 2,  37, 30, 14, 4,  31, 24, 15, 8,  38,
                            26, 17, 6,  33, 23, 13, 5,  34, 29, 18,
                            7,  36, 28, 20, 3,  35, 21, 16, 1,  32,
                            25, 12, 10, 40, 22, 11, 9,  39, 27, 19 },
                          "random ranks match the pinned values");
}

void
test_fenwick_tree()
{
  for (size_t n : { 0, 1, 2, 7, 64, 100 }) {
    wdm::utils::FenwickTree tree(n);
    std::vector<double> values(n, 0.0);
    bool agrees = tree.prefix_sum(0) == 0.0;
    for (size_t step = 0; step < 3 * n; ++step) {
      size_t index = (step * 37) % n;
      double value = (step % 5 == 0) ? -1.5 : 0.5 * static_cast<double>(step);
      tree.add(index, value);
      values[index] += value;
      for (size_t end = 0; end <= n; ++end) {
        double expected = 0.0;
        for (size_t k = 0; k < end; ++k)
          expected += values[k];
        agrees = agrees && std::fabs(tree.prefix_sum(end) - expected) < 1e-9;
      }
    }
    test::check(agrees, "FenwickTree prefix sums match the running sums");
  }
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

void
test_soft_ranks()
{
  const std::vector<std::string> methods{ "min", "average", "first", "random" };
  const double scale = 1.5e-8;

  // without near-ties, the soft ranks are the ranks by value, bit for bit
  std::vector<double> x{ 3, 1, 2, 2, 5, 3, 2, 4, 1, 2 };
  std::vector<double> w{ 1, 2, 1, 1, 3, 1, 2, 1, 1, 1 };
  bool same = true;
  for (const auto& method : methods) {
    for (const auto& weights : { std::vector<double>(), w }) {
      same = same && (wdm::impl::rank(x, weights, method, { 7 }, scale) ==
                      wdm::impl::rank(x, weights, method, { 7 }));
    }
  }
  test::check(same, "soft ranks equal the ranks by value without near-ties");

  // a crowd of distinct values within rounding of each other: shifting it by
  // far less than the scale moves the soft ranks by a negligible amount, where
  // the ranks by value reorder it
  std::vector<double> a, b;
  for (int i = 0; i < 400; ++i) {
    const double v = static_cast<double>(i % 50) / 50.0;
    a.push_back(v + static_cast<double>((i % 7) - 3) * 1e-13);
    b.push_back(v + static_cast<double>((i % 5) - 2) * 1e-13);
  }
  for (const auto& method : { "average", "first", "random" }) {
    const auto ra = wdm::impl::rank(a, {}, method, { 3 }, scale);
    const auto rb = wdm::impl::rank(b, {}, method, { 3 }, scale);
    const auto ha = wdm::impl::rank(a, {}, method, { 3 });
    const auto hb = wdm::impl::rank(b, {}, method, { 3 });
    double soft_move = 0.0, hard_move = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      soft_move = std::max(soft_move, std::fabs(ra[i] - rb[i]));
      hard_move = std::max(hard_move, std::fabs(ha[i] - hb[i]));
    }
    test::check(soft_move < 1e-3,
                std::string("soft ranks move continuously, ") + method);
    test::check(hard_move >= 1.0,
                std::string("ranks by value reorder the crowd, ") + method);
    // every pair splits its weight between the two orders
    double total = 0.0;
    for (double r : ra)
      total += r;
    test::check_near(total,
                     400.0 * 401.0 / 2.0,
                     std::string("soft ranks sum to n(n + 1) / 2, ") + method,
                     1e-8);
  }

  // missing values stay missing, and the scale must be a distance
  auto with_nan = wdm::impl::rank(
    { 0.5, std::numeric_limits<double>::quiet_NaN(), 0.5 + 1e-12 },
    {},
    "average",
    {},
    scale);
  test::check(std::isnan(with_nan[1]) && !std::isnan(with_nan[0]) &&
                (std::fabs(with_nan[0] + with_nan[2] - 3.0) < 1e-12),
              "soft ranks skip missing values");
  for (double bad : { -1.0,
                      std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::quiet_NaN() })
    test::check_throws([&]() { wdm::impl::rank(x, {}, "average", {}, bad); },
                       "rank rejects a scale that is not a finite distance");
}
}

int
main()
{
  test_rank_ties();
  test_order_within_ties();
  test_tie_keys();
  test_random_ties_are_local();
  test_random_ranks_are_portable();
  test_rank0_ties();
  test_soft_ranks();
  test_fenwick_tree();
  return test::finish();
}
