/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "csc.hh"
#include "pagerankpipeline.hh"
#include <cstdio>
#include <cmath>

int main () {
  // Run an example on kernel 3.  This example was computed using the
  // Octave version on 4/1/2016.
  typedef std::tuple<uint32_t, uint32_t, double> triple;
  std::vector<triple> tuples = {
    triple(1-1, 2-1, 0.30000),
    triple(2-1, 2-1, 0.25000),
    triple(3-1, 2-1, 0.23077),
    triple(5-1, 2-1, 0.18182),
    triple(7-1, 2-1, 0.50000),
    triple(8-1, 2-1, 1),
    triple(1-1, 3-1, 0.26667),
    triple(2-1, 3-1, 0.16667),
    triple(3-1, 3-1, 0.15385),
    triple(5-1, 3-1, 0.27273),
    triple(6-1, 3-1, 1),
    triple(1-1, 4-1, 0.066667),
    triple(7-1, 4-1, 0.50000),
    triple(1-1, 5-1, 0.20000),
    triple(2-1, 5-1, 0.25000),
    triple(3-1, 5-1, 0.23077),
    triple(4-1, 5-1, 0.50000),
    triple(5-1, 5-1, 0.27273),
    triple(1-1, 6-1, 0.066667),
    triple(3-1, 6-1, 0.23077),
    triple(5-1, 6-1, 0.090909),
    triple(1-1, 7-1, 0.10000),
    triple(2-1, 7-1, 0.16667),
    triple(3-1, 7-1, 0.15385),
    triple(5-1, 7-1, 0.18182),
    triple(2-1, 8-1, 0.16667),
    triple(4-1, 8-1, 0.50000)
  };
  csc_matrix<uint32_t> nonzeros(8, tuples);
  std::vector<double> initial_r = {0.076423, 0.485566, 0.454543,  0.794698,  0.047933, 0.653576, 0.681263, 0.934423};
  std::vector<double> r = kernel3_compute<uint32_t>(3, nonzeros, &initial_r);
  std::vector<double> expected_r = {0.018750, 0.334375, 0.191262, 0.089931, 0.152989, 0.042293, 0.053178, 0.117222};
  if (0) {
    printf("r = {");
    for (auto v : r) {
      printf("%f, ", v);
    }
    printf("}\n");
  }
  assert(expected_r.size() == r.size());
  //printf("diff = {");
  for (size_t i = 0; i < r.size(); i++) {
    double ri = r[i];
    double ei = expected_r[i];
    double relerr = std::abs(ri-ei) / std::max(std::abs(ri), std::abs(ei));
    //printf(" %g (%g),\n", r[i]-expected_r[i], relerr);
    assert(relerr < 1e-5);
  }
  //printf("}\n");
  return 0;
}

// Local Variables:
// mode: C++
// c-basic-offset: 2
// indent-tabs-mode: nil
// End:
