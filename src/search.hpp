#pragma once

#include "geometry.hpp"

#include <chrono>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace polyjam {

struct SearchConfig {
  int count = 1;
  double seconds = 10.0;
  int threads = 1;
  std::uint64_t seed = 1;
  double start_scale = 0.0;
  double tolerance = 1e-5;
  double clearance = 0.0;
};

struct Metrics {
  double violation = 0.0;
  double max_violation = 0.0;
  double containment_sq = 0.0;
  double overlap_sq = 0.0;
  double hull_volume = 0.0;
  int overlapping_pairs = 0;
};

struct Packing {
  std::vector<Pose> poses;
  double scale = 0.0;
  Metrics metrics;
  std::uint64_t iterations = 0;
  double elapsed = 0.0;
  bool feasible = false;
};

using ProgressCallback = std::function<void(const Packing&)>;

Metrics evaluate(const Polyhedron& piece, const Polyhedron& shell, const Packing& p, double clearance = 0.0);
Packing search(const Polyhedron& piece, const Polyhedron& shell, const SearchConfig& config, ProgressCallback progress = {});
void write_result_json(
  const std::string& path,
  const Polyhedron& piece,
  const Polyhedron& shell,
  const SearchConfig& config,
  const Packing& result,
  const std::string& note = {}
);
std::string progress_json(const Packing& p, const char* type);

} // namespace polyjam
