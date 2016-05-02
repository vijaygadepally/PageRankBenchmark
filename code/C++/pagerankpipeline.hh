/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#ifndef PAGERANKPIPELINE_HH
#define PAGERANKPIPELINE_HH

#include "csc.hh"
#include <cstdio>
#include <vector>

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files);

extern FILE *data_file; // set to non-NULL if you want a data file suitable for gnuplot.

// Internal functions exported for testing:
template <class T>
std::vector<double> kernel3_compute(const int SCALE, 
                                    const csc_matrix<T> &M,
                                    // for testing we use a known r.
                                    std::vector<double> *initial_r = nullptr);

#endif
