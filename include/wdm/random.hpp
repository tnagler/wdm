
#pragma once

#include <algorithm> // For std::generate
#include <cstdint>
#include <random> // For std::mt19937, std::seed_seq, std::random_device
#include <vector>

namespace wdm {

namespace random {

//! Random-number generator used for reproducible randomized tie breaking.
//!
//! The same draws on every platform: the engine (`std::mt19937`) and its
//! seeding (`std::seed_seq`) are specified exactly by the standard, and the
//! distributions are implemented here rather than taken from the standard
//! library, whose distributions are implementation-defined.
class RandomGenerator
{
public:
  //! @param seeds seeds of the generator; if empty, it is seeded randomly.
  explicit RandomGenerator(std::vector<int> seeds = std::vector<int>())
    : generator(initialize_generator(seeds))
  {
  }

  //! draws a size_t uniformly in [0, n - 1]; `n` must be positive.
  size_t sample_int(size_t n)
  {
    if (static_cast<uint64_t>(n) > UINT32_MAX) {
      return static_cast<size_t>(sample_wide(static_cast<uint64_t>(n)));
    }
    // Lemire's multiply-shift: the high half of draw * range is uniform once
    // draws whose low half falls below 2^32 mod range are rejected
    const uint32_t range = static_cast<uint32_t>(n);
    uint64_t product = static_cast<uint64_t>(generator()) * range;
    uint32_t low = static_cast<uint32_t>(product);
    if (low < range) {
      const uint32_t threshold = (0u - range) % range;
      while (low < threshold) {
        product = static_cast<uint64_t>(generator()) * range;
        low = static_cast<uint32_t>(product);
      }
    }
    return static_cast<size_t>(product >> 32);
  }

  //! draws a double uniformly in [0, 1), on a grid of 2^-53.
  double sample_double()
  {
    return static_cast<double>(next64() >> 11) / 9007199254740992.0;
  }

private:
  std::mt19937 generator;

  uint64_t next64()
  {
    const uint64_t high = static_cast<uint64_t>(generator());
    return (high << 32) | static_cast<uint64_t>(generator());
  }

  // `sample_int` past 2^32 - 1: 2^64 mod range is rejected, which leaves a
  // multiple of `range` values, so the remainder is exactly uniform
  uint64_t sample_wide(uint64_t range)
  {
    const uint64_t threshold = (0 - range) % range;
    uint64_t draw;
    do {
      draw = next64();
    } while (draw < threshold);
    return draw % range;
  }

  static std::mt19937 initialize_generator(std::vector<int>& seeds)
  {
    if (seeds.empty()) {
      seeds = generate_random_seeds();
    }
    std::seed_seq seq(seeds.begin(), seeds.end());
    return std::mt19937(seq);
  }

  // Generate random seeds using std::random_device
  static std::vector<int> generate_random_seeds()
  {
    std::random_device rd{};
    std::vector<int> seeds(5);
    std::generate(
      seeds.begin(), seeds.end(), [&]() { return static_cast<int>(rd()); });
    return seeds;
  }
};

// Custom shuffle function
template<typename T>
void
shuffle(std::vector<T>& vec, RandomGenerator& rand_gen)
{
  for (size_t i = vec.size() - 1; i > 0; --i) {
    size_t j = rand_gen.sample_int(i + 1); // Generate a random index in [0, i]
    std::swap(vec[i], vec[j]);
  }
}

}
}
