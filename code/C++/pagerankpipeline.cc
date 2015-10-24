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
void read_files(int kernel, int SCALE, int edges_per_vertex, int n_files, 
                std::vector<std::tuple<T, T>> *edges) {
  const uint64_t M_per_file = (1u<<SCALE) * edges_per_vertex;
  const uint64_t M = M_per_file * n_files;
  edges->reserve(M);
  for (int i = 0; i < n_files; i++) {
    FILE *f = fopen(file_named(kernel, SCALE, i).c_str(), "r");
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
      edges->push_back(std::tuple<T,T>(v1, v2));
    }
    if (line) free(line);
    fclose(f);
  }
  assert(edges->size() == M);
}


template <class T>
void write_files(int kernel, int SCALE, int edges_per_vertex, int n_files, 
                 const std::vector<std::tuple<T, T>> &edges) {
  const uint64_t M_per_file = (1u<<SCALE) * edges_per_vertex;
  auto it = edges.begin();
  mkdir(dir_named(1, SCALE).c_str(), 0777);
  for (int i = 0; i < n_files; i++) {
    std::ofstream f(file_named(kernel, SCALE, i), std::ios::out);
    for (T j = 0; j < M_per_file; j++) {
      f << std::get<0>(*it) << '\t' << std::get<1>(*it) << '\n';
      it++;
    }
  }
  assert(it == edges.end());

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
void kernel1(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  // Sort the data
  std::vector<std::tuple<T, T>> edges;
  read_files<T>(0, SCALE, edges_per_vertex, n_files, &edges);
  std::sort(edges.begin(), edges.end());
  write_files<T>(1, SCALE, edges_per_vertex, n_files, edges);
  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K1time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * edges.size() / (end-start));
}

template <class T>
void kernel2(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  std::vector<std::tuple<T, T>> edges;
  read_files<T>(1, SCALE, edges_per_vertex, n_files, &edges);

  // Remove duplicates and construct a matrix.
  // Construct the transposed matrix to make it easy to remove the supernode and leaves.
  std::vector<std::tuple<T, T, double>> matrix;
  matrix.reserve(edges.size());
  {
    auto &prev = edges[0];
    double count = 1;
    for (size_t i = 1; i < edges.size(); i++) {
      if (prev == edges[i]) {
        count++;
      } else {
        // store column then row (the matrix is transposed)
        matrix.push_back(std::tuple<T, T, double>(std::get<1>(prev), std::get<0>(prev), count));
        prev = edges[i];
        count = 1;
      }
    }
    // store column then row (the matri is transposed)
    matrix.push_back(std::tuple<T, T, double>(std::get<1>(prev), std::get<0>(prev), count));
  }

  // Remove the supernodes and leaves.  The in-degree is the number of nonzeros in a column.
  // Get rid of the column(s) with the most entries  and the columns with exactly one entry.
  // Since we stored columns we can sort it.
  std::sort(matrix.begin(), matrix.end());

  // Find the maximum column
  size_t max_col_count = std::numeric_limits<size_t>::max();
  {
    T      prev_col = std::get<0>(matrix[0]);
    size_t prev = 0;
    for (size_t i = 1; i < matrix.size(); i++) {
      // if it's the same column we keep going
      if (std::get<0>(matrix[i]) == prev_col) {
        continue;
      } else {
        // it's a new column
        size_t n_in_col = i - prev;
        if (n_in_col > max_col_count) {
          max_col_count = n_in_col;
        }
        prev = i;
      }
    }
    size_t n_in_col = matrix.size() - prev;
    if (n_in_col > max_col_count) {
      max_col_count = n_in_col;
    }
  }
  
  // Untranpose the matrix
  {
    size_t write_here = 0;
    T      prev_col = std::get<0>(matrix[0]);
    size_t prev = 0;
    for (size_t i = 1; i < matrix.size(); i++) {
      // If it's the same column we keep going
      if (std::get<0>(matrix[i]) == prev_col) {
        continue;
      } else {
        // it's a new column
        size_t n_in_col = i - prev;
        if (n_in_col != 1 && n_in_col != max_col_count) {
          // copy the column to the new matrix.
          while (prev < i) {
            const auto &e = matrix[prev];
            matrix[write_here++] = std::tuple<T,T,double>(std::get<1>(e), std::get<0>(e), std::get<2>(e));
            prev++;
          }
          prev = i;
        }
      }
    }
    size_t n_in_col = matrix.size() - prev;
    if (n_in_col != 1 && n_in_col != max_col_count) {
      // copy the column to the new matrix.
      while (prev < matrix.size()) {
        const auto &e = matrix[prev];
        matrix[write_here++] = std::tuple<T,T,double>(std::get<1>(e), std::get<0>(e), std::get<2>(e));
        prev++;
      }
    }
    matrix.resize(write_here);
  }


  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K2time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * edges.size() / (end-start));
}

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  kernel0<T>(SCALE, edges_per_vertex, n_files);
  kernel1<T>(SCALE, edges_per_vertex, n_files);
  kernel2<T>(SCALE, edges_per_vertex, n_files);
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
