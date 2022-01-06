# wdm

![build status](https://github.com/tnagler/wdm/actions/workflows/main.yml/badge.svg?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> A header-only C++ library for weighted dependence measures

Provides efficient implementations of weighted dependence measures and related 
independence tests:

- Pearsons's rho
- Spearmans's rho
- Kendall's tau
- Blomqvist's beta
- Hoeffding's D

All measures are computed in O(_n log n_) time, where _n_ is the number of 
observations.

### Functionality

The library provides:

- a function `wdm()` to compute the weighted dependence measures,
- a class `Indep_test` to perform a test for independence based on asymptotic
  p-values.

For details, see the [API documentation](https://tnagler.github.io/wdm/) 
and the [example](#example) below.

### Dependencies

The library only requires C++11. 

For projects already using the [Eigen](https://eigen.tuxfamily.org) linear 
algebra library, there are convenience wrappers that can be made available via 

```cpp
#include <wdm/eigen.hpp>
```

### Including the library in other projects

There are two options: 

1. Either copy the header files in `include/` to your project.
2. Install the headers globally using the CMake project. To do that go to the 
   root repository of this repo and run:
   ```shell
   mkdir build && cd build         # open build folder
   cmake .. && sudo make install   # install library
   cd .. && rm -rf build           # leave and remove build folder
   ```
   To use the library in your project, just add 
   `target_link_libraries(your_proj_name wdm)` to `your_proj_name/CMakeLists.txt`.

You can then include the main header in your source code:

```cpp
#include <wdm.hpp>
```

### Example

```cpp
#include "wdm.hpp"

// input vectors
std::vector<double> x{1, 3, 2, 5, 3, 2, 20, 15};
std::vector<double> y{2, 12, 4, 7, 8, 14, 17, 6};

// weights
std::vector<double> w{1, 1, 2, 2, 1, 0, 0.5, 0.3};

std::cout <<
    "unweighted Kendall's tau: " << wdm::wdm(x, y, "kendall") << std::endl;
std::cout <<
    "weighted Kendall's tau: " <<  wdm::wdm(x, y, "kendall", w) << std::endl;

// weighted independence test
wdm::Indep_test test(x, y, "kendall", w);
std::cout << "statistic: " << test.statistic() << std::endl;
std::cout << "p-value: " << test.p_value() << std::endl;
```

```
unweighted Kendall's tau: 0.2965
weighted Kendall's tau: 0.550633
statistic: 1.71047
p-value: 0.0871793
```
