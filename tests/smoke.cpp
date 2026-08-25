#include "geometry.hpp"

#include <cmath>
#include <iostream>

using namespace polyjam;

static void require(bool ok, const char* message) {
  if (!ok) { std::cerr << "FAIL: " << message << "\n"; std::exit(1); }
}

int main() {
  const auto cube = builtin_polyhedron("cube");
  const auto dodeca = builtin_polyhedron("dodeca");
  require(cube.vertices.size() == 8, "cube vertex count");
  require(cube.faces.size() == 6, "cube face count");
  require(cube.edges.size() == 12, "cube edge count");
  require(dodeca.vertices.size() == 20, "dodeca vertex count");
  require(dodeca.faces.size() == 12, "dodeca face count");

  Pose origin{};
  const auto v0 = transformed_vertices(cube, origin);
  const auto contained = containment_violation(cube, 1.0, v0);
  require(contained.max_excess < 1e-10, "identity cube contained at scale 1");
  const auto too_small = containment_violation(cube, 0.9, v0);
  require(too_small.max_excess > 0.01, "scale .9 must violate cube containment");

  Pose a{};
  Pose b{};
  a.position = {-0.7,0,0};
  b.position = { 0.7,0,0};
  const auto av = transformed_vertices(cube,a);
  const auto bv = transformed_vertices(cube,b);
  const auto separated = sat_test(cube,a,av,b,bv);
  require(separated.separated, "separated cubes detected");

  b.position = {0.3,0,0};
  const auto bv2 = transformed_vertices(cube,b);
  const auto overlapping = sat_test(cube,a,av,b,bv2);
  require(!overlapping.separated, "overlapping cubes detected");
  require(overlapping.penetration > 0.01, "overlap penetration positive");

  std::cout << "geometry smoke tests passed\n";
}
