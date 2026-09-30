
#pragma once

#include <algorithm> // For std::generate
#include <cstdint>
#include <random> // For std::mt19937, std::seed_seq, std::random_device
#include <vector>

namespace wdm {

namespace random {

//! Random-number generator used for reproducible randomized tie breaking.
//!
//! The same draws on every platform, with or without Boost: the engine
//! (`std::mt19937`) and its seeding (`std::seed_seq`) are specified exactly by
//! the standard, and the distributions are implemented here rather than taken
//! from the standard library, whose distributions are implementation-defined.
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
    const uint64_t range = static_cast<uint64_t>(n);
    // 2^64 mod range: accepting only draws at or above it leaves a multiple
    // of `range` values, so the remainder is exactly uniform
    const uint64_t threshold = (0 - range) % range;
    uint64_t draw;
    do {
      draw = next64();
    } while (draw < threshold);
    return static_cast<size_t>(draw % range);
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
