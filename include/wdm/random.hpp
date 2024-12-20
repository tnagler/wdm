
#include <vector>

#ifdef USE_BOOST
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/seed_seq.hpp>
#else
#include <algorithm> // For std::generate
#include <random>    // For std::random_device
#endif

namespace wdm
{

  namespace random
  {

    class RandomGenerator
    {
    public:
      // Constructor with optional seeds
      explicit RandomGenerator(std::vector<int> seeds = std::vector<int>())
#ifdef USE_BOOST
          : generator(initialize_boost_generator(seeds)){}
#else
          : generator(initialize_std_generator(seeds))
      {
      }
#endif

            // Sample a size_t in [0, n-1]
            size_t sample_int(size_t n)
      {
#ifdef USE_BOOST
        boost::random::uniform_int_distribution<size_t> distribution(0, n - 1);
#else
        std::uniform_int_distribution<size_t> distribution(0, n - 1);
#endif
        return distribution(generator);
      }

      // Sample a double in [0.0, 1.0)
      double sample_double()
      {
#ifdef USE_BOOST
        boost::random::uniform_real_distribution<double> distribution(0.0, 1.0);
#else
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
#endif
        return distribution(generator);
      }

    private:
#ifdef USE_BOOST
      boost::random::mt19937 generator;

      // Initialize Boost generator with seeds
      boost::random::mt19937 initialize_boost_generator(std::vector<int> &seeds)
      {
        if (seeds.empty())
        {
          seeds = generate_random_seeds();
        }
        boost::random::seed_seq seq(seeds.begin(), seeds.end());
        return boost::random::mt19937(seq);
      }
#else
      std::default_random_engine generator;

      // Initialize std generator with seeds
      std::default_random_engine initialize_std_generator(std::vector<int> &seeds)
      {
        if (seeds.empty())
        {
          seeds = generate_random_seeds();
        }
        std::seed_seq seq(seeds.begin(), seeds.end());
        return std::default_random_engine(seq);
      }
#endif

      // Generate random seeds using std::random_device
      static std::vector<int> generate_random_seeds()
      {
        std::random_device rd{};
        std::vector<int> seeds(5);
        std::generate(seeds.begin(), seeds.end(), [&]()
                      { return static_cast<int>(rd()); });
        return seeds;
      }
    };

    // Custom shuffle function
    template <typename T>
    void shuffle(std::vector<T> &vec, RandomGenerator &rand_gen)
    {
      for (size_t i = vec.size() - 1; i > 0; --i)
      {
        size_t j = rand_gen.sample_int(i + 1); // Generate a random index in [0, i]
        std::swap(vec[i], vec[j]);
      }
    }

  }
}