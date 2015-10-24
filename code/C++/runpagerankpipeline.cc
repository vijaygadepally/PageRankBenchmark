/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#include "pagerankpipeline.hh"

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <limits>


static const int min_scale_default = 10;
static const int max_scale_default = 16;

static char *progname;

static void usage() __attribute__((noreturn));
static void usage() {
  printf("Usage: %s [ min_scale [ max_scale ]]\n", progname);
  printf("  Runs the page rank pipeline on a range of scales\n");
  printf("  low_scale defaults to %d\n", min_scale_default);
  printf("  hi_scale defaults to %d\n", max_scale_default);
  exit(1);
}

static unsigned int parse_uint(const char *str) {
  errno = 0;
  char *end;
  long long v = strtol(str, &end, 10);
  if (errno != 0 || end == str || *end != 0
      || v < 0
      || v > std::numeric_limits<unsigned int>::max()) {
    printf("Unable to parse this as an unsigned integer: %s", str);
    usage();
  }
  return v;
}

int main(int argc, char *argv[] __attribute__((unused))) {
  progname = argv[0];
  const int min_scale = (argc >= 2) ? parse_uint(argv[1]) : min_scale_default;
  const int max_scale = (argc >= 3) ? parse_uint(argv[2]) : max_scale_default;
  if (argc >= 4) {
    printf("Too many arguments\n");
    usage();
  }
  const int edges_per_vertex = 6;
  const int nfile = 4;
  for (int scale = min_scale; scale <= max_scale; scale++) {
    pagerankpipeline<uint32_t>(scale, edges_per_vertex, nfile);
  }
  return 0;
}
