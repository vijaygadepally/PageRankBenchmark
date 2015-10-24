/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "pagerankpipeline.hh"

#include "fasttime.h"
#include "krongraph500.hh"

#include <iostream>
#include <fstream>

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  for (int i = 0; i < n_files; i++) {
    std::ofstream f("data/K" + std::to_string(SCALE) + "/" + std::to_string(i) + ".tsv", std::ios::out);
    for (const auto &pair : kronecker<T>(SCALE, edges_per_vertex, false)) {
      f << std::get<0>(pair) << '\t' << std::get<1>(pair) << '\n';
    }
    // The stream closes itself when leaving scope.
  }
  fasttime_t end   = gettime();
  printf("K0time: %.6fs, edges/sec: %g\n", end-start, ((1u<<SCALE)*edges_per_vertex*n_files)/(end-start));
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
