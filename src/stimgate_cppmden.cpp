/*
 * stimgate_cppmden.cpp — cpp11 R interface for the FAUST-derived cpPmden
 * taut-string density estimator.
 *
 * Responsibilities:
 *   1. Accept a numeric vector from R.
 *   2. Remove non-finite values at the R level (done before calling this).
 *   3. Sort the values (cpPmden requires sorted input).
 *   4. Call cpPmden().
 *   5. Return the `stringInfo.string` density vector and `nmax`.
 *
 * This file is part of stimgate and is licensed under GPL (>= 3).
 * The underlying cpPmden() implementation is derived from FAUST
 * (faster annotation using shape-constrained trees), which is also
 * GPL (>= 3).  See inst/COPYRIGHTS for full provenance details.
 */

#include <cpp11.hpp>
#include <vector>
#include <algorithm>
#include "faust.h"

[[cpp11::register]]
cpp11::list stimgate_cpPmden(std::vector<double> x) {
  // Sort: cpPmden() assumes sorted input.
  std::sort(x.begin(), x.end());

  unsigned long n = x.size();

  // Return degenerate result for inputs too small for the algorithm.
  if (n < 3) {
    cpp11::writable::doubles dens(n > 1 ? n - 1 : 0);
    for (auto& v : dens) v = 0.0;
    using cpp11::literals::operator""_nm;
    return cpp11::writable::list({"string"_nm = dens, "nmax"_nm = 0});
  }

  stringInfo res = cpPmden(x);

  cpp11::writable::doubles dens(res.string.begin(), res.string.end());
  using cpp11::literals::operator""_nm;
  return cpp11::writable::list({"string"_nm = dens, "nmax"_nm = res.nmax});
}
