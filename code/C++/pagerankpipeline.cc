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
void write_files(const int kernel, const int SCALE, const int edges_per_vertex, const int n_files, 
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
void kernel0(const int SCALE, const int edges_per_vertex, const int n_files) {
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
void kernel1(const int SCALE, const int edges_per_vertex, const int n_files) {
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
struct csr_matrix {
  std::vector<T> row_starts;
  std::vector<T> cols;
  std::vector<double> vals;
};

template <class T>
void kernel2(const int SCALE, const int edges_per_vertex, const int n_files, csr_matrix<T> *result) {
  fasttime_t start = gettime();
  const size_t N = (1u<<SCALE);
  std::vector<std::tuple<T, T>> edges;
  read_files<T>(1, SCALE, edges_per_vertex, n_files, &edges);

  // Remove duplicates and construct a matrix.
  // Construct the transposed matrix to make it easy to remove the supernode and leaves.
  std::vector<std::tuple<T, T, double>> matrix;
  matrix.reserve(edges.size());
  {
    auto &prev = edges[0];
    T count = 1;
    for (size_t i = 1; i < edges.size(); i++) {
      if (prev == edges[i]) {
        count++;
      } else {
        // store column then row (the matrix is transposed)
        matrix.push_back(std::tuple<T, T, T>(std::get<1>(prev), std::get<0>(prev), count));
        prev = edges[i];
        count = 1;
      }
    }
    // store column then row (the matri is transposed)
    matrix.push_back(std::tuple<T, T, T>(std::get<1>(prev), std::get<0>(prev), count));
  }

  // Remove the supernodes and leaves.  The in-degree is the number of nonzeros in a column.
  // Get rid of the column(s) with the most entries  and the columns with exactly one entry.
  // Since we stored columns we can sort it.
  std::sort(matrix.begin(), matrix.end());

  {
    std::vector<T> col_counts(N, 0);
    for (const auto &e : matrix) {
      col_counts[ std::get<0>(e) ] += std::get<2>(e);
    }
    if (0) {
      printf("col counts =");
      for (const auto &c : col_counts) {
        printf(" %d", c);
      }
      printf("\n");
    }
    T max_col_count = *std::max_element(col_counts.begin(), col_counts.end());
    //printf("Max col count = %d\n", max_col_count);

    std::vector<double> col_inverses(N, 0);
    for (size_t i = 0; i < N; i++) {
      if (col_counts[i] != 0) col_inverses[i] = 1/(double)(col_counts[i]);
    }

    // Untranpose the matrix and remove the supernode and leaves.
    size_t new_size = 0;
    for (const auto &e : matrix) {
      const T &col = std::get<0>(e);
      const T &row = std::get<1>(e);
      const double &val = std::get<2>(e);
      const T col_count = col_counts[col];
      if (col_count == 1 || col_count == max_col_count) continue;
      matrix[new_size++] = std::tuple<T,T,double>(row, col, val*col_inverses[col]);
    }
    matrix.resize(new_size);
    matrix.shrink_to_fit();
  }

  // Put it in row-major order.
  std::sort(matrix.begin(), matrix.end());
  if (0) {
    for (const auto &e : matrix) {
      printf(" (%d %d %f)", std::get<0>(e), std::get<1>(e), std::get<2>(e));
    }
    printf("\n");
  }

  // Construct the csr matrix.
  auto &row_starts = result->row_starts;
  auto &cols       = result->cols;
  auto &vals       = result->vals;
  row_starts.reserve(N + 1);
  cols.reserve(matrix.size());
  vals.reserve(matrix.size());
  for (const auto & e : matrix) {
    const T &row = std::get<0>(e);
    const T &col = std::get<1>(e);
    const double &val = std::get<2>(e);
    while (row_starts.size() <= static_cast<size_t>(row)) {
      row_starts.push_back(cols.size());
    }
    cols.push_back(col);
    vals.push_back(val);
  }
  while (row_starts.size() < N+1) {
    row_starts.push_back(cols.size());
  }
  if (0)
  for (size_t row = 0; row < N; row++) {
    double sum = 0;
    for (size_t i = row_starts[row]; i < row_starts[row+1]; i++) {
      sum += vals[i];
    }
    double inv = 1/sum;
    for (size_t i = row_starts[row]; i < row_starts[row+1]; i++) {
      vals[i] *= inv;
    }
  }
  if (0) {
    printf("row_starts:");
    for (auto &rs : row_starts) { printf(" %d", rs); }
    printf("\ncols: ");
    for (auto &cs : cols) { printf(" %4d", cs); }
    printf("\nvals: ");
    for (auto &vs : vals) { printf(" %3.3f", vs); }
    printf("\n");
  }

  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K2time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * edges.size() / (end-start));
}

template <class T>
void kernel3(const int SCALE, const int edges_per_vertex, const csr_matrix<T> &M __attribute__((unused))) {
  fasttime_t start = gettime();
  const size_t N = (1u<<SCALE);
  std::vector<double> r;
  r.reserve(N);
  uint64_t sum = 0;
  for (size_t i = 0; i < N; i++) {
    long rval = random();
    r.push_back((double)rval);  
    sum += rval;
  }
  double one_over_sum = 1/(double)sum;
  double fsum=0;
  for (auto &e : r) {
    e *= one_over_sum;
    fsum += e;
  }
  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K2time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * (1u<<SCALE)*edges_per_vertex / (end-start));
}

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  kernel0<T>(SCALE, edges_per_vertex, n_files);
  kernel1<T>(SCALE, edges_per_vertex, n_files);
  csr_matrix<T> M;
  kernel2<T>(SCALE, edges_per_vertex, n_files, &M);
  kernel3<T>(SCALE, edges_per_vertex,           M);
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
