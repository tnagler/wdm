# wdm
A header-only C++11 library for weighted dependence measures

Provides an efficient implementation of weighted dependence measures:

   * Pearsons's rho: `prho()`
   * Spearmans's rho: `srho()`
   * Kendall's tau: `ktau()`
   * Blomqvist's beta: `bbeta()`
   * Hoeffding's D: `hoeffd()`

All functions take at least two arguments for the two input vectors to compute 
the dependence measure on, plus an optional argument to supply weights for 
each observation. Additionally asymptotic p-values for all measures can be 
calculated with `indep_test()`.

### Example

``` cpp
#include "wdm.hpp"

// input vectors
std::vector<double> x{1, 3, 2, 5, 3, 2, 20, 15};
std::vector<double> y{2, 12, 4, 7, 8, 14, 17, 6};

// weights
std::vector<double> w{1, 1, 2, 2, 1, 0, 0.5, 0.3};

// unweighted Kendall's tau
std::cout << "Kendall's tau: " << wdm::ktau(x, y) << std::endl;
// or: wdm::wdm(x, y, "ktau")
std::cout << "p-value: " << wdm::indeptest(x, y, "ktau") << std::endl;

// weighted Kendall's tau
std::cout << "Kendall's tau: " << wdm::ktau(x, y, w) << std::endl;
// or: wdm::wdm(x, y, "ktau", w)
std::cout << "p-value: " << wdm::indeptest(x, y, "ktau", w) << std::endl;
```
Output:
```
Kendall's tau: 0.2965
p-value: 0.208413
Kendall's tau: 0.550633
p-value: 0.0557333
```

### How to use the library

There are two options: 

1. Either copy the header files in your project.
2. Install the headers globally using the CMake project. To do that go to the 
   root repository of this repo and run:
   ```shell
   mkdir build && cd build         # open build folder
   cmake .. && sudo make install   # install library
   cd .. && rm -rf build           # leave and remove build folder
   ```

You can then include the header containing all functionality
``` cpp
#include <wdm.hpp>
``` 
or include only specific headers, e.g.
``` cpp
#include <wdm/ktau.hpp>
``` 

### Dependencies

The library only requires C++11 and the STL. For projects already using the [Eigen](https://eigen.tuxfamily.org) linear algebra library, there are 
convenience wrappers that can be made available via 
``` cpp
#include <wdm/eigen.hpp>
``` 
