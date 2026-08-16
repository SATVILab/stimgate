#ifndef STIMGATE_FAUST_H
#define STIMGATE_FAUST_H

/*
 * Minimal header for the FAUST-derived taut-string / cpPmden implementation
 * vendored into stimgate.
 *
 * The original faust.h is part of FAUST (faster annotation using
 * shape-constrained trees), Copyright (C) RGLab, licensed under
 * GPL (>= 3).  Only the declarations actually required by the vendored
 * source files are reproduced here; all unrelated FAUST machinery
 * (annotation/clustering structs, tsGates, kMedDP, singleDip, etc.)
 * has been omitted.
 *
 * See inst/COPYRIGHTS for full provenance information.
 */

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>

// Result structure returned by tautString() and cpPmden().
// `string` is the piecewise-constant density (length n-1 for n input points).
struct stringInfo {
  std::vector<double> string;
  std::vector<int>    knotsind;
  std::vector<double> knotst;
  std::vector<double> knotsy;
  int nknots;
  int nmax;
};

// Forward declarations of the functions used across the vendored files.
void local_density(const std::vector<double>&, std::vector<int>&,
                   long, long, long);

stringInfo tautString(const std::vector<double>&,
                      const std::vector<double>&,
                      const std::vector<double>&,
                      const std::vector<double>&,
                      double, double, long, int);

void easymax(std::vector<double>&, long, long, long,
             long*, long*, double*);

void difficultmax(std::vector<double>&, long, long, long,
                  long*, long*, double*);

std::vector<double> kkuiper(std::vector<double>&, long, int);

std::vector<double> cppApprox(std::vector<double>&,
                               std::vector<double>&,
                               std::vector<double>&);

stringInfo cpPmden(const std::vector<double>&);

std::vector<double> rQuantile(const std::vector<double>&,
                               std::vector<double>);

double medianAbsoluteDeviation(const std::vector<double>&);

#endif // STIMGATE_FAUST_H
