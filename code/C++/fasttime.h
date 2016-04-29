/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/* Pagerank Pipeline Benchmark in C++                          */
/* Copyright 2015 Bradley C. Kuszmaul, bradley@mit.edu         */

#ifndef FASTTIME_H
#define FASTTIME_H

#include <assert.h>

#ifdef __MACH__
#include <mach/mach_time.h> // mach_absolute_time

typedef uint64_t fasttime_t;

static inline fasttime_t gettime(void)
// Effect: return the current time.
{
  return mach_absolute_time();
}

static inline double tdiff(fasttime_t start, fasttime_t end)
// Effect: Return the time different between the start and the end, as a float
//  in units of seconds.  This function does not need to be fast.
// Implementation notes:  See https://developer.apple.com/library/mac/qa/qa1398/_index.html
{
  static mach_timebase_info_data_t timebase;
  int r = mach_timebase_info(&timebase);
  assert(r == 0);
  fasttime_t elapsed = end-start;
  double ns = (double)elapsed * timebase.numer / timebase.denom;
  return ns*1e-9;
}

static inline unsigned int random_seed_from_clock(void) {
  fasttime_t now = gettime();
  return (now & 0xFFFFFFFF) + (now>>32);
}

#else // LINUX

#include <time.h>

typedef struct timespec fasttime_t;

static inline fasttime_t gettime(void)
// Effect: return the current time.
{
  struct timespec s;
  int r = clock_gettime(CLOCK_MONOTONIC, &s);
  assert(r == 0);
  return s;
}

static inline double tdiff(fasttime_t start, fasttime_t end)
// Effect: Return the time different between the start and the end, as a float
//  in units of seconds.  This function does not need to be fast.
{
  return end.tv_sec - start.tv_sec + 1e-9*(end.tv_nsec - start.tv_nsec);
}

static inline unsigned int random_seed_from_clock(void) {
  fasttime_t now = gettime();
  return now.tv_sec + now.tv_nsec;
}

// Poison these symbols to help find portability problems.
int clock_gettime(clockid_t, struct timespec *) __attribute__((deprecated));
time_t time(time_t *) __attribute__((deprecated));

#endif // LINUX

#ifdef __cplusplus
static inline double operator- (const fasttime_t &end, const fasttime_t &start) {
    return tdiff(start, end);
}
#endif

#endif //FASTTIME_H
