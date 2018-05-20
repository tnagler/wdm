# wdm
A header-only C++11 library for weighted dependence measures

Provides efficient implementations of weighted dependence measures:

   * Pearsons's rho
   * Spearmans's rho
   * Kendall's tau
   * Blomqvist's beta
   * Hoeffding's D

### Usage

The core functions, `wdm()` and `indep_test()`, take the following arguments:

   * `x`, `y`: the input vectors to compute the dependence measure on.
   * `method`: A `std::string` indicating the method.
   * `weights`: a vector of weights for each observation (optional).

`wdm()` computes the dependence measure, and `indep_test()` calculates the 
asymptotic p-value under the Null hypothesis of independence. For more 
details, see the [API documentation](https://tnagler.github.io/wdm/).


### Example

``` cpp
#include "wdm.hpp"

// input vectors
std::vector<double> x{1, 3, 2, 5, 3, 2, 20, 15};
std::vector<double> y{2, 12, 4, 7, 8, 14, 17, 6};

// weights
std::vector<double> w{1, 1, 2, 2, 1, 0, 0.5, 0.3};

// unweighted Kendall's tau
std::cout << "Kendall's tau: " << wdm::wdm(x, y, "ktau") << std::endl;
// or: wdm::ktau(x, y)
std::cout << "p-value: " << wdm::indep_test(x, y, "ktau") << std::endl;

// weighted Kendall's tau
std::cout << "Kendall's tau: " <<  wdm::wdm(x, y, "ktau", w) << std::endl;
// or: wdm::ktau(x, y, w)
std::cout << "p-value: " << wdm::indep_test(x, y, "ktau", w) << std::endl;
```
Output:
```
Kendall's tau: 0.2965
p-value: 0.208413
Kendall's tau: 0.550633
p-value: 0.0557333
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

You can then include the header containing all functionality
``` cpp
#include <wdm.hpp>
``` 
or include only specific headers, e.g.
``` cpp
#include <wdm/ktau.hpp>
``` 

### Dependencies

The library only requires C++11 and [Boost](https://www.boost.org/). 

For projects already using the [Eigen](https://eigen.tuxfamily.org) linear 
algebra library, there are convenience wrappers that can be made available via 
``` cpp
#include <wdm/eigen.hpp>
``` 
See the [API documentation](https://tnagler.github.io/wdm/).
