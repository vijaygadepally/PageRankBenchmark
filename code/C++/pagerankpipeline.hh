/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#ifndef PAGERANKPIPELINE_HH
#define PAGERANKPIPELINE_HH

#include <cstdio>

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files);

extern FILE *data_file; // set to non-NULL if you want a data file suitable for gnuplot.


#endif
