#include "geometry.hpp"
#include "search.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace polyjam;

int main() {
  const auto cube = builtin_polyhedron("cube");
  Packing single_cube;
  single_cube.scale = 1.0;
  single_cube.poses.resize(1);
  const Metrics single_metrics = evaluate(cube, cube, single_cube);
  const double expected_cube_volume = 8.0 / (3.0 * std::sqrt(3.0));
  if (std::abs(single_metrics.hull_volume - expected_cube_volume) > 1e-8) {
    std::cerr << "FAIL: configuration hull volume for one cube\n";
    return EXIT_FAILURE;
  }
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

  const auto dodeca = builtin_polyhedron("dodeca");
  const Packing dodeca_result = search(cube, dodeca, cfg);
  if (!dodeca_result.feasible || dodeca_result.scale > 4.2 || dodeca_result.poses.size() != 50) {
    std::cerr << "FAIL: expected a dense feasible 50-cube dodecahedron construction\n";
    return EXIT_FAILURE;
  }

  const auto tetra = builtin_polyhedron("tetra");
  cfg.count = 30;
  const Packing tetra_result = search(tetra, dodeca, cfg);
  if (!tetra_result.feasible || tetra_result.scale > 3.5 || tetra_result.poses.size() != 30) {
    std::cerr << "FAIL: expected a dense feasible tetrahedron lattice construction\n";
    return EXIT_FAILURE;
  }

  const auto icosa = builtin_polyhedron("icosa");
  cfg.count = 17;
  cfg.seed = 19;
  const Packing threshold_result = search(tetra, icosa, cfg);
  if (!threshold_result.feasible || threshold_result.scale > 3.1 || threshold_result.poses.size() != 17) {
    std::cerr << "FAIL: expected a feasible tetrahedron crystal-core threshold construction\n";
    return EXIT_FAILURE;
  }

  cfg.count = 50;
  cfg.seed = 7;
  const Packing round_result = search(icosa, dodeca, cfg);
  if (!round_result.feasible || round_result.scale > 4.8 || round_result.poses.size() != 50) {
    std::cerr << "FAIL: expected a dense sphere-like lattice construction for icosahedra\n";
    return EXIT_FAILURE;
  }
  std::cout << "lattice search test passed\n";
}
