/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "pagerankpipeline.hh"

#include "fasttime.h"
#include "krongraph500.hh"
#include "csc.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#include <sys/stat.h>

FILE *data_file = NULL; // set to non-null to produce a data file suitable for gnuplot.



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
  const uint64_t M_per_file = (1ul<<SCALE) * edges_per_vertex;
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
  const uint64_t M_per_file = (1ul<<SCALE) * edges_per_vertex;
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
         1e-6 * ((1ul<<SCALE)*edges_per_vertex*n_files) / (end-start),
         bytes_written,
         1e-6 * bytes_written / (end-start));
  if (data_file) {
    fprintf(data_file, "%d %ld %g", SCALE, (1ul<<SCALE)*edges_per_vertex, ((1ul<<SCALE)*edges_per_vertex*n_files) / (end-start));
  }
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
  if (data_file) {
    fprintf(data_file, " %g", edges.size() / (end-start));
  }
}

template <class T>
csc_matrix<T> kernel2(const int SCALE, const int edges_per_vertex, const int n_files) {
  fasttime_t start = gettime();
  const size_t N = (1ul<<SCALE);
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

    // Untranspose the matrix and remove the supernode and leaves.
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

  // Don't put it into row-major order, since we want to put the matrix into CSC format, column-major sort order is good.

  // Construct the csr matrix.
  csc_matrix<T> result(N, matrix);
  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K2time: %9.3fs Medges/sec: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * edges.size() / (end-start));
  if (data_file) {
    fprintf(data_file, " %g", edges.size() / (end-start));
  }
  return result;
}

const int page_rank_iteration_count = 1;

// Fast random number generator.  Don't need high-quality.
// A hard-coded pseudorandom number generator
static uint64_t X = 0xce3f12500545b241ul;
static uint64_t a = 0x3b812326f4736033ul;
static uint64_t b = 0x0a58ec3022b3b942ul;

inline uint64_t
prandnum() {
  X = a*X + b;
  return X;
}

template <class T>
std::vector<double> kernel3_compute(const int SCALE, 
                                    const csc_matrix<T> &M,
                                    // for testing we use a known r.
                                    std::vector<double> *initial_r) {

  const size_t N = (1ul<<SCALE);
  std::vector<double> r;
  double sum = 0; // this will compute norm(r, 1), also written in matab as sum(r,2)

  if (initial_r != nullptr) {
    assert(initial_r->size() == N);
    r = *initial_r;
    for (double rval : r) sum += rval;
  } else {
    // Generate a random starting rank.  This corresponds to matlab's
    //       r = rand(1, N)
    r.reserve(N);
    for (size_t i = 0; i < N; i++) {
      double rval = (prandnum()>>32);
      r.push_back(rval);
      sum += rval;
    }
  }
  if (0) {
    printf("r = ");
    for (const auto &v : r) {
      printf("%f ", v);
    }
    printf("\n");
    printf("sum=%f\n", sum);
  }

  // Normalize.  This corresponds to matlab's r = r ./ norm(r,1)
  double one_over_sum = 1/(double)sum;
  double fsum=0;
  for (auto &e : r) {
    e *= one_over_sum;
    fsum += e;
  }

  if (0) {
    printf("r = ");
    for (auto & e : r) {
      printf(" %f", e);
    }
    printf("\n");
  }

  // Now do the page rank.
  double c = 0.85;

  // In matlab, a is a vector of these values, but we need represent
  //   the value of a only once.
  double a = (1-c)/(double)N;

  // In matlab we write r = ((c .* r) * M) + (a .* sum(r,2))
  // In C++, we create a second vector r2, and std::swap r2 with r at
  //  the end.
  std::vector<double> r2(N, 0);
  for (int pr_count = 0; pr_count < page_rank_iteration_count; pr_count++) {
    for (size_t i = 0; i < N; i++) {
      // In matlab, this is    r = ((c .* r) * M) + (a .* sum(r,2))
      double dotsum = 0;
      const T start_col = M.col_starts[i];
      const T end_col   = M.col_starts[i+1];
      for (T vi = start_col; vi < end_col; vi++) {
        dotsum += r[M.rows[vi]] * M.vals[vi];
      }
      r2[i] = c * dotsum  + a * fsum;

#if 0
      // In matlab, this is:   r =  M * (R .* c) + a;
      double sum = 0;
      //printf("Row %ld\n", i);
      const T start_row = M.row_starts[i];
      const T end_row   = M.row_starts[i+1];
      for (T vi = start_row; vi < end_row; vi++) {
        //printf("sum += %f * %f\n", M.vals[vi], r[M.cols[vi]]);
        sum += M.vals[vi] * r[M.cols[vi]];
      }
      r2[i] = sum * c + a;
#endif
    }
    std::swap(r, r2);
    if (0) {
      printf("after iteration %d r=", pr_count);
      for (auto & e : r) {
        printf(" %f", e);
      }
      printf("\n");
    }
  }
  return r;
}

template <class T>
void kernel3(const int SCALE, const int edges_per_vertex, const csc_matrix<T> &M) {
  fasttime_t start = gettime();
  std::vector<double> r = kernel3_compute<T>(SCALE, M);
  fasttime_t end   = gettime();
  printf("scale=%2d Edgefactor=%2d K3time: %9.3fs Medges/sec: %7.2f  MFLOPS: %7.2f\n", 
         SCALE, edges_per_vertex, 
         end-start,
         1e-6 * (1ul<<SCALE)*edges_per_vertex*page_rank_iteration_count / (end-start),
         2e-6 * (1ul<<SCALE)*edges_per_vertex*page_rank_iteration_count / (end-start));
  if (data_file) {
    fprintf(data_file, " %g\n", (1ul<<SCALE)*edges_per_vertex*page_rank_iteration_count / (end-start));
  }
}

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  kernel0<T>(SCALE, edges_per_vertex, n_files);
  kernel1<T>(SCALE, edges_per_vertex, n_files);
  csc_matrix<T> M = kernel2<T>(SCALE, edges_per_vertex, n_files);
  kernel3<T>(SCALE, edges_per_vertex,           M);
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
template std::vector<double> kernel3_compute<uint32_t>(int, csc_matrix<uint32_t> const&, std::vector<double> *);

// Local Variables:
// mode: C++
// c-basic-offset: 2
// indent-tabs-mode: nil
// End:
