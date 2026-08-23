/*
 * stimgate_cppmden.cpp — cpp11 R interface for the FAUST-derived cpPmden
 * taut-string density estimator.
 *
 * Responsibilities:
 *   1. Accept a numeric vector from R.
 *   2. Remove non-finite values defensively before entering native code.
 *   3. Sort the values (cpPmden requires sorted input).
 *   4. Return a benign zero-density result for degenerate inputs.
 *   5. Call cpPmden().
 *   6. Return the `stringInfo.string` density vector and `nmax`.
 *
 * This file is part of stimgate and is licensed under GPL (>= 3).
 * The underlying cpPmden() implementation is derived from FAUST
 * (faster annotation using shape-constrained trees), which is also
 * GPL (>= 3).  See inst/COPYRIGHTS for full provenance details.
 */

#include <cpp11.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include "faust.h"

namespace {
cpp11::list zero_density_result(unsigned long n) {
  cpp11::writable::doubles dens(n > 1 ? n - 1 : 0);
  std::fill(dens.begin(), dens.end(), 0.0);
  using cpp11::literals::operator""_nm;
  return cpp11::writable::list({"string"_nm = dens, "nmax"_nm = 0});
}
}  // namespace

[[cpp11::register]]
cpp11::list stimgate_cpPmden(std::vector<double> x) {
  x.erase(
    std::remove_if(
      x.begin(),
      x.end(),
      [](double value) { return !std::isfinite(value); }
    ),
    x.end()
  );

  // Sort: cpPmden() assumes sorted input.
  std::sort(x.begin(), x.end());

  unsigned long n = x.size();

  // cpPmden uses a Kuiper-bound table starting at n = 50; smaller inputs
  // may not converge. Constant input cannot define a density over positive-
  // width intervals and would otherwise create zero denominators downstream.
  if (n < 50) {
    return zero_density_result(n);
  }

  double range = x.back() - x.front();
  if (!std::isfinite(range) || range <= 0.0) {
    return zero_density_result(n);
  }

  stringInfo res = cpPmden(x);

  cpp11::writable::doubles dens(res.string.begin(), res.string.end());
  using cpp11::literals::operator""_nm;
  return cpp11::writable::list({"string"_nm = dens, "nmax"_nm = res.nmax});
}
