/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "pagerankpipeline.hh"

#include "fasttime.h"
#include "krongraph500.hh"

#include <iostream>
#include <fstream>

#include <sys/stat.h>

static std::string dir_named(int SCALE)  {
  return "data/K" + std::to_string(SCALE);
}
  

static std::string file_named(int SCALE, int filenum) {
  return dir_named(SCALE) + "/" + std::to_string(filenum) + ".tsv";
}

template <class T>
void pagerankpipeline(int SCALE, int edges_per_vertex, int n_files) {
  fasttime_t start = gettime();
  mkdir("data", 0777);
  for (int i = 0; i < n_files; i++) {
    mkdir(dir_named(SCALE).c_str(), 0777);
    if (0) {
      std::ofstream f(file_named(SCALE, i), std::ios::out);
      for (const auto &pair : kronecker<T>(SCALE, edges_per_vertex, false)) {
        f << std::get<0>(pair) << '\t' << std::get<1>(pair) << '\n';
      }
      // The stream closes itself when leaving scope.
    } else {
      write_kronecker(file_named(SCALE, i), SCALE, edges_per_vertex);
    }
  }
  fasttime_t end   = gettime();
  uint64_t bytes_written = 0;
  for (int i = 0; i < n_files; i++) {
    struct stat statbuf;
    int r = stat(file_named(SCALE, i).c_str(), &statbuf);
    assert(r == 0);
    bytes_written += statbuf.st_size;
  }
  printf("Scale=%d Edgefactor=%d K0time: %9.3fs, Medges/sec: %6.2f  bytes:%12ld Mbytes/sec: %6.2f\n", 
         SCALE, edges_per_vertex,
         end-start,
         1e-6 * ((1u<<SCALE)*edges_per_vertex*n_files) / (end-start),
         bytes_written,
         1e-6 * bytes_written / (end-start));
}

template void pagerankpipeline<uint32_t>(int SCALE, int edges_per_vertex, int n_files);
