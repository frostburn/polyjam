#include "geometry.hpp"
#include "search.hpp"

#include <cstdlib>
#include <iostream>

using namespace polyjam;

int main() {
  const auto cube = builtin_polyhedron("cube");
  SearchConfig cfg;
  cfg.count = 50;
  cfg.seconds = 0.02;
  cfg.threads = 1;
  cfg.seed = 7;

  const Packing result = search(cube, cube, cfg);
  if (!result.feasible || result.scale > 4.0000001 || result.poses.size() != 50) {
    std::cerr << "FAIL: expected a feasible 50-cube construction at scale 4\n";
    return EXIT_FAILURE;
  }
  std::cout << "cube lattice search test passed\n";
}
