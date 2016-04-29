/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#ifndef KRONGRAPH500_HH
#define KRONGRAPH500_HH

#include <tuple>
#include <vector>

template <class T>
std::vector<std::tuple<T, T>> kronecker(int SCALE, int edges_per_vertex, bool randomize);
// Effect: Create an edge list according to the Graph500 randomized Kronecker graph.
//   SCALE is the log (base 2) of the total number of vertice.s
//   edges_per_vertex is the average number of edges per node
//   randomize: if true, randomize the permutation, otherwise don't.
//   Duplicate edges may be returned, all the dges ha
//  Requires that ((1<<SCALE) * edges_per_vertex) fits in T.

void write_kronecker(const std::string &filename, int SCALE, int edges_per_verte);
// Effect: Write an edgelist to a file.  For this we don't have the randomize option (we never randomize).
//  Edges are written as tab-separated integers followed by a newline.

#endif
