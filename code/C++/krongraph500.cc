/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "krongraph500.hh"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>

#include <fcntl.h>
#include <unistd.h>

template <class T>
std::vector<std::tuple<T, T>> kronecker(int SCALE, int edges_per_vertex, bool randomize) {
  // Half-hearted attempt to check that the T is big enough.
  // Doesn't try do a good job if SCALE and edges_per_vertex are close to 1ul<<64.
  assert(std::numeric_limits<T>::max() >= (1ul << SCALE)*edges_per_vertex);
  T N = T(1)<<SCALE;
  T M = edges_per_vertex * N;
  double A = 0.57, B = 0.19, C = 0.19;
  double ab = A+B;
  double c_norm = C/(1 - (A + B));
  double a_norm = A/(A+B);
  std::vector<std::tuple<T,T>> edges;
  edges.reserve(M);
  uint64_t rand_max_scaled = RAND_MAX;
  uint64_t ab_scaled = rand_max_scaled * ab;
  uint64_t c_norm_scaled = rand_max_scaled * c_norm;
  uint64_t a_norm_scaled = rand_max_scaled * a_norm;
  for (T i = 0; i < M; i++) {
    T ij_i = 0, ij_j = 0;
    for (int ib = 0; ib < SCALE; ib++) {
      T r1 = random();
      T r2 = random();
      T ii_bit = r1 > ab_scaled;
      T jj_bit = r2 > (c_norm_scaled * ii_bit + a_norm_scaled * !ii_bit);
      ij_i += ii_bit << ib;
      ij_j += jj_bit << ib;
    }
    edges.push_back(std::tuple<T,T>(ij_i, ij_j));
  }
  if (randomize) {
    // Permute the vertex labels.
    {
      std::vector<T> p(N);
      std::iota(p.begin(), p.end(), 0);
      std::random_shuffle(p.begin(), p.end());
      for (auto &e : edges) {
        e = std::tuple<T, T>(p[std::get<0>(e)],
                             p[std::get<1>(e)]);
      }
    }
    // Permute the edge list
    {
      std::random_shuffle(edges.begin(), edges.end());
    }
  }
  return edges;
}

template std::vector<std::tuple<uint32_t, uint32_t>> kronecker(int SCALE, int edges_per_vertex, bool randomize);

static int appendint_internal(int bufoff, char *buf, uint64_t v) {
  if (v<10) {
    buf[bufoff++] = '0'+v;
    return bufoff;
  } else {
    return appendint_internal(appendint_internal(bufoff, buf, v/10),
                              buf,
                              v%10);
  }
}
static int appendint(int bufoff, char *buf, uint64_t v) {
  if (v==0) {
    buf[bufoff++] = '0';
    return bufoff;
  } else {
    return appendint_internal(bufoff, buf, v);
  }
}

void write_kronecker(const std::string &filename, int SCALE, int edges_per_vertex) {
  // Half-hearted attempt to check that the T is big enough.
  // Doesn't try do a good job if SCALE and edges_per_vertex are close to 1ul<<64.
  int fd = open(filename.c_str(), O_WRONLY | O_TRUNC | O_CREAT, 0777);
  assert(fd >= 0);
  assert(std::numeric_limits<uint64_t>::max() >= (1ul << SCALE)*edges_per_vertex);
  uint64_t N = uint64_t(1)<<SCALE;
  uint64_t M = edges_per_vertex * N;
  double A = 0.57, B = 0.19, C = 0.19;
  double ab = A+B;
  double c_norm = C/(1 - (A + B));
  double a_norm = A/(A+B);
  uint64_t rand_max_scaled = RAND_MAX;
  uint64_t ab_scaled = rand_max_scaled * ab;
  uint64_t c_norm_scaled = rand_max_scaled * c_norm;
  uint64_t a_norm_scaled = rand_max_scaled * a_norm;
  const int bufsize = 1024 * 1024;
  const int buf_full_mark = bufsize-1024;
  char *buf = new char [bufsize];
  int bufoff = 0;
  for (uint64_t i = 0; i < M; i++) {
    uint64_t ij_i = 0, ij_j = 0;
    for (int ib = 0; ib < SCALE; ib++) {
      uint64_t r1 = random();
      uint64_t r2 = random();
      uint64_t ii_bit = r1 > ab_scaled;
      uint64_t jj_bit = r2 > (c_norm_scaled * ii_bit + a_norm_scaled * !ii_bit);
      ij_i += ii_bit << ib;
      ij_j += jj_bit << ib;
    }
    bufoff = appendint(bufoff, buf, ij_i);
    buf[bufoff++] = '\t';
    bufoff = appendint(bufoff, buf, ij_j);    
    buf[bufoff++] = '\n';
    if (bufoff > buf_full_mark) {
      ssize_t ws = write(fd, buf, bufoff);
      assert(ws == bufoff);
      bufoff = 0;
    }
  }
  if (bufoff > 0) {
    ssize_t ws = write(fd, buf, bufoff);
    assert(ws == bufoff);
    bufoff = 0;
  }
  close(fd);
  delete [] buf;
}


