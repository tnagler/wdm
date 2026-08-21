# wdm

![build status](https://github.com/tnagler/wdm/actions/workflows/main.yml/badge.svg?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`wdm` is a header-only C++11 library implementing weighted dependence
measures and related asymptotic independence tests. It primarily provides the
computational core for higher-level interfaces; users of the R interface
should consult that interface's documentation for end-user workflows.

All estimators have an average time complexity of O(_n_ log _n_).

## Supported methods

| Method | Accepted names | Direction | Minimum sample |
| --- | --- | --- | ---: |
| Pearson correlation | `pearson`, `prho`, `cor` | symmetric | 2 |
| Spearman's rho | `spearman`, `srho`, `rho` | symmetric | 2 |
| Kendall's tau | `kendall`, `ktau`, `tau` | symmetric | 2 |
| Blomqvist's beta | `blomqvist`, `bbeta`, `beta` | symmetric | 2 |
| Hoeffding's D | `hoeffding`, `hoeffd`, `d` | symmetric | 5 |
| Chatterjee's xi | `chatterjee`, `cxi`, `xi` | `y` on `x` | 2 |

The primary C++ entry points are:

- `wdm::wdm()` for an estimate;
- `wdm::Indep_test` for an estimate, test statistic, effective sample size,
  and asymptotic p-value;
- overloads in `<wdm/eigen.hpp>` for Eigen vectors and matrices.

See the [API documentation](https://tnagler.github.io/wdm/) for signatures.

## Input behavior

Input vectors must have equal sizes. Optional case weights must be finite,
nonnegative, and have a positive total. Multiplying all weights by a positive
constant does not change an estimate or its inference, and zero-weight rows
are ignored.

By default, rows containing `NaN` in either variable or the weights are
removed. If too few observations remain, estimates and inference results are
`NaN`. With `remove_missing = false`, missing or insufficient input raises
`std::runtime_error`. Size mismatches, invalid weights, and unknown method or
alternative names also raise `std::runtime_error`.

`Indep_test` supports `two-sided`, `less`, and `greater` alternatives, except
that Hoeffding's D is two-sided only. Its weighted approximations use Kish's
effective sample size and require enough effective observations for the
selected transformation.

### Chatterjee's xi

Chatterjee's xi measures dependence of `y` on `x`, so reversing its arguments
can change the result. Ties in `x` are broken uniformly at random without
consulting `y`; pass `seeds` to reproduce the same tie ordering.

The weighted estimate supports tied responses. Analytic inference with unequal
weights currently requires a continuous response, weights that are fixed or
depend only on `x`, and sufficiently diffuse normalized weights. Set
`y_continuous` to `false` when the response distribution is discrete, even if
the observed sample has no ties. Observed response ties select the discrete
inference path automatically. Unequally weighted inference is unavailable in
either discrete case.

For a continuous response, `estimate()` returns the general
denominator-corrected coefficient while `statistic()` standardizes the
continuous-response inferential approximation.

## Using the C++ headers

No compiled library is required. The only mandatory dependency is C++11.
Either copy `include/` into a project or install the CMake package:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build
cmake --install build --prefix /desired/prefix
```

Consumers can then use:

```cmake
find_package(wdm CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE wdm)
```

Set `CMAKE_PREFIX_PATH` when installing to a nonstandard prefix. Include the
main API with:

```cpp
#include <wdm.hpp>
```

Eigen convenience overloads require Eigen and are enabled by including
`<wdm/eigen.hpp>`. The standard-library random backend is the default; configure
with `-DUSE_BOOST=ON` to use Boost.Random instead.

## Example

```cpp
#include <iostream>
#include <wdm.hpp>

int main()
{
  std::vector<double> x{ 1, 3, 2, 5, 3, 2, 20, 15 };
  std::vector<double> y{ 2, 12, 4, 7, 8, 14, 17, 6 };
  std::vector<double> weights{ 1, 1, 2, 2, 1, 0, 0.5, 0.3 };

  std::cout << "unweighted Kendall's tau: "
            << wdm::wdm(x, y, "kendall") << '\n';
  std::cout << "weighted Kendall's tau: "
            << wdm::wdm(x, y, "kendall", weights) << '\n';

  wdm::Indep_test test(x, y, "kendall", weights);
  std::cout << "statistic: " << test.statistic() << '\n';
  std::cout << "p-value: " << test.p_value() << '\n';
}
```

The example prints, to the shown precision:

```text
unweighted Kendall's tau: 0.2965
weighted Kendall's tau: 0.550633
statistic: 1.41025
p-value: 0.158465
```

## Development

The project uses the style in `.clang-format`. Check tracked C++ files with:

```sh
git ls-files -z '*.h' '*.hpp' '*.ipp' '*.cpp' '*.cc' |
  xargs -0 clang-format --dry-run --Werror
```

Regenerate the checked-in API documentation from the repository root with:

```sh
doxygen docs/Doxyfile.in
```
