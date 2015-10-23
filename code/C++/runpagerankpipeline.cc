#include "pagerankpipeline.hh"

#include <cassert>
#include <cstdint>

int main(int argc, char *argv[] __attribute__((unused))) {
    assert(argc == 1);
    const int scale_min = 10;
    const int scale_max = 10;
    const int edges_per_vertex = 6;
    const int nfile = 4;
    for (int scale = scale_min; scale <= scale_max; scale++) {
        pagerankpipeline<uint32_t>(scale, edges_per_vertex, nfile);
    }
    return 0;
}
