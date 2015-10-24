/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "pagerankpipeline.hh"

#include "fasttime.h"
#include "krongraph500.hh"

#include <algorithm>
#include <iostream>
#include <fstream>

#include <sys/stat.h>

static std::string dir_named(int kernel, int SCALE)  {
  return "data/K" + std::to_string(kernel) + "_S" + std::to_string(SCALE);
}
  

static std::string file_named(int kernel, int SCALE, int filenum) {
  return dir_named(kernel, SCALE) + "/" + std::to_string(filenum) + ".tsv";
}


template <class T>
void kernel0(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  mkdir("data", 0777);
  for (int i = 0; i < n_files; i++) {
    mkdir(dir_named(0, SCALE).c_str(), 0777);
    if (0) {
      std::ofstream f(file_named(0, SCALE, i), std::ios::out);
      for (const auto &pair : kronecker<T>(SCALE, edges_per_vertex, false)) {
        f << std::get<0>(pair) << '\t' << std::get<1>(pair) << '\n';
      }
      // The stream closes itself when leaving scope.
    } else {
      write_kronecker(file_named(0, SCALE, i), SCALE, edges_per_vertex);
    }
  }
  fasttime_t end   = gettime();
  uint64_t bytes_written = 0;
  for (int i = 0; i < n_files; i++) {
    struct stat statbuf;
    int r = stat(file_named(0, SCALE, i).c_str(), &statbuf);
    assert(r == 0);
    bytes_written += statbuf.st_size;
  }
  printf("Scale=%2d Edgefactor=%2d K0time: %9.3fs Medges/sec: %7.2f  bytes:%12ld Mbytes/sec: %6.2f\n", 
         SCALE, edges_per_vertex,
         end-start,
         1e-6 * ((1u<<SCALE)*edges_per_vertex*n_files) / (end-start),
         bytes_written,
         1e-6 * bytes_written / (end-start));
}

template <class T>
T parse_int(char *str, char **end) {
  uint64_t r = 0;
  while (1) {
    char c = *str;
    if (c >= '0' && c <= '9') {
      r = r*10 + c - '0';
      str++;
    } else {
      *end = str;
      return r;
    }
  }
}
template <class T>
void kernel1(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  // Sort the data
  std::vector<std::tuple<T, T>> edges;
  const uint64_t M_per_file = (1u<<SCALE) * edges_per_vertex;
  const uint64_t M = M_per_file * n_files;
  edges.reserve(M);
  for (int i = 0; i < n_files; i++) {
    FILE *f = fopen(file_named(0, SCALE, i).c_str(), "r");
    assert(f);
    char *line = NULL;
    size_t len = 0;
    ssize_t n_read;
    while ((n_read = getline(&line, &len, f)) != -1) {
      char *tab;
      T v1 = parse_int<T>(line, &tab);
      assert(tab[0] == '\t');
      char *nl;
      T v2 = parse_int<T>(tab+1, &nl);
      assert(nl[0] == '\n');
      assert(nl[1] == 0);
      edges.push_back(std::tuple<T,T>(v1, v2));
    }
    if (line) free(line);
    fclose(f);
  }
  assert(edges.size() == M);
  std::sort(edges.begin(), edges.end());
  auto it = edges.begin();
  mkdir(dir_named(1, SCALE).c_str(), 0777);
  for (int i = 0; i < n_files; i++) {
    std::ofstream f(file_named(1, SCALE, i), std::ios::out);
    for (T j = 0; j < M_per_file; j++) {
      f << std::get<0>(*it) << '\t' << std::get<1>(*it) << '\n';
      it++;
    }
  }
  assert(it == edges.end());
  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K1time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * M / (end-start));
}
template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  kernel0<T>(SCALE, edges_per_vertex, n_files);
  kernel1<T>(SCALE, edges_per_vertex, n_files);
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
