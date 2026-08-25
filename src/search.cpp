#include "search.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace polyjam {

namespace {

struct HullTriangle { int a; int b; int c; };

double configuration_hull_volume(const std::vector<std::vector<Vec3>>& pieces) {
  std::vector<Vec3> points;
  for (const auto& vertices : pieces) points.insert(points.end(), vertices.begin(), vertices.end());
  if (points.size() < 4) return 0.0;

  const int p0 = static_cast<int>(std::min_element(points.begin(), points.end(), [](const Vec3& a, const Vec3& b) {
    return a.x < b.x;
  }) - points.begin());
  int p1 = -1;
  double best = 0.0;
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    const double d = norm2(points[static_cast<std::size_t>(i)] - points[static_cast<std::size_t>(p0)]);
    if (d > best) { best = d; p1 = i; }
  }
  if (p1 < 0 || best < 1e-20) return 0.0;

  int p2 = -1;
  best = 0.0;
  const Vec3 line = points[static_cast<std::size_t>(p1)] - points[static_cast<std::size_t>(p0)];
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    const double d = norm2(cross(line, points[static_cast<std::size_t>(i)] - points[static_cast<std::size_t>(p0)]));
    if (d > best) { best = d; p2 = i; }
  }
  if (p2 < 0 || best < 1e-20) return 0.0;

  int p3 = -1;
  best = 0.0;
  const Vec3 plane = cross(points[static_cast<std::size_t>(p1)] - points[static_cast<std::size_t>(p0)],
                           points[static_cast<std::size_t>(p2)] - points[static_cast<std::size_t>(p0)]);
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    const double d = std::abs(dot(plane, points[static_cast<std::size_t>(i)] - points[static_cast<std::size_t>(p0)]));
    if (d > best) { best = d; p3 = i; }
  }
  if (p3 < 0 || best < 1e-14) return 0.0;

  const Vec3 interior = (points[static_cast<std::size_t>(p0)] + points[static_cast<std::size_t>(p1)] +
                         points[static_cast<std::size_t>(p2)] + points[static_cast<std::size_t>(p3)]) / 4.0;
  auto outward = [&](HullTriangle f) {
    const Vec3& a = points[static_cast<std::size_t>(f.a)];
    const Vec3 n = cross(points[static_cast<std::size_t>(f.b)] - a, points[static_cast<std::size_t>(f.c)] - a);
    if (dot(n, interior - a) > 0.0) std::swap(f.b, f.c);
    return f;
  };
  std::vector<HullTriangle> faces{
    outward({p0,p1,p2}), outward({p0,p3,p1}), outward({p0,p2,p3}), outward({p1,p3,p2})
  };

  for (int p = 0; p < static_cast<int>(points.size()); ++p) {
    std::vector<bool> visible(faces.size(), false);
    bool any = false;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      const auto& f = faces[i];
      const Vec3& a = points[static_cast<std::size_t>(f.a)];
      const Vec3 n = cross(points[static_cast<std::size_t>(f.b)] - a, points[static_cast<std::size_t>(f.c)] - a);
      visible[i] = dot(n, points[static_cast<std::size_t>(p)] - a) > 1e-10 * std::max(1.0, norm(n));
      any = any || visible[i];
    }
    if (!any) continue;
    std::vector<std::pair<int,int>> horizon;
    auto add_edge = [&](int a, int b) {
      const auto reverse = std::find(horizon.begin(), horizon.end(), std::pair<int,int>{b,a});
      if (reverse != horizon.end()) horizon.erase(reverse); else horizon.emplace_back(a,b);
    };
    std::vector<HullTriangle> kept;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      if (!visible[i]) kept.push_back(faces[i]);
      else { add_edge(faces[i].a, faces[i].b); add_edge(faces[i].b, faces[i].c); add_edge(faces[i].c, faces[i].a); }
    }
    faces = std::move(kept);
    for (const auto& edge : horizon) faces.push_back(outward({edge.first, edge.second, p}));
  }

  double volume = 0.0;
  for (const auto& f : faces) volume += dot(points[static_cast<std::size_t>(f.a)],
    cross(points[static_cast<std::size_t>(f.b)], points[static_cast<std::size_t>(f.c)])) / 6.0;
  return std::abs(volume);
}

} // namespace

Metrics evaluate(const Polyhedron& piece, const Polyhedron& shell, const Packing& p, double clearance) {
  Metrics m;
  const int n = static_cast<int>(p.poses.size());
  std::vector<std::vector<Vec3>> world;
  world.reserve(static_cast<std::size_t>(n));
  for (const auto& pose : p.poses) world.push_back(transformed_vertices(piece, pose));

  for (int i=0;i<n;i++) {
    const auto c = containment_violation(shell, p.scale, world[static_cast<std::size_t>(i)]);
    m.containment_sq += c.sum_sq;
    m.max_violation = std::max(m.max_violation, c.max_excess);
  }

  for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
    const auto sep = sat_test(piece, p.poses[static_cast<std::size_t>(i)], world[static_cast<std::size_t>(i)],
                             p.poses[static_cast<std::size_t>(j)], world[static_cast<std::size_t>(j)]);
    if (!sep.separated) {
      const double pen = sep.penetration + clearance;
      if (pen > 0.0) {
        m.overlap_sq += pen*pen;
        m.max_violation = std::max(m.max_violation, pen);
        ++m.overlapping_pairs;
      }
    }
  }
  m.violation = m.containment_sq + m.overlap_sq;
  m.hull_volume = configuration_hull_volume(world);
  return m;
}

static double rand01(std::mt19937_64& rng) {
  return std::generate_canonical<double, 53>(rng);
}
static double normal01(std::mt19937_64& rng) {
  static thread_local std::normal_distribution<double> dist(0.0,1.0);
  return dist(rng);
}
static Vec3 random_unit(std::mt19937_64& rng) {
  Vec3 v{normal01(rng),normal01(rng),normal01(rng)};
  if (norm2(v) < 1e-20) v = {1,0,0};
  return normalized(v);
}
static Quat random_q(std::mt19937_64& rng) {
  const double u1=rand01(rng), u2=rand01(rng), u3=rand01(rng);
  const double a=std::sqrt(1.0-u1), b=std::sqrt(u1);
  const double t1=2.0*std::numbers::pi*u2, t2=2.0*std::numbers::pi*u3;
  return normalized({b*std::cos(t2),a*std::sin(t1),a*std::cos(t1),b*std::sin(t2)});
}

static Vec3 random_center_in_shell(const Polyhedron& shell, double scale, std::mt19937_64& rng) {
  // Hull is normalized to circumradius 1, so [-scale,scale]^3 contains it.
  for (int attempt=0; attempt<5000; ++attempt) {
    Vec3 p{(2*rand01(rng)-1)*scale, (2*rand01(rng)-1)*scale, (2*rand01(rng)-1)*scale};
    bool inside = true;
    for (const auto& f : shell.faces) if (dot(f.normal,p) > scale*f.d) { inside=false; break; }
    if (inside) return p;
  }
  return {};
}

static Packing random_initial(const Polyhedron& shell, const SearchConfig& cfg, double scale, std::mt19937_64& rng) {
  Packing p;
  p.scale = scale;
  p.poses.resize(static_cast<std::size_t>(cfg.count));
  for (auto& pose : p.poses) {
    pose.position = random_center_in_shell(shell, scale*0.85, rng);
    pose.rotation = random_q(rng);
  }
  return p;
}

static Packing cube_grid_initial(const SearchConfig& cfg, std::mt19937_64& rng) {
  // A random cloud is an especially poor starting point for cube-in-cube
  // instances: at useful densities almost every cube overlaps several others.
  // Start those searches from an exact lattice construction instead.  Search
  // all integer box dimensions because, for non-perfect cubes, a rectangular
  // grid can be substantially tighter than ceil(cuberoot(n)) in every axis.
  int best_x = cfg.count;
  int best_y = 1;
  int best_z = 1;
  int best_side = cfg.count;
  int best_capacity = cfg.count;
  for (int x = 1; x <= cfg.count; ++x) {
    for (int y = x; y <= cfg.count; ++y) {
      const int xy = x * y;
      const int z = (cfg.count + xy - 1) / xy;
      if (z < y) continue;
      const int side = z;
      const int capacity = xy * z;
      if (side < best_side || (side == best_side && capacity < best_capacity)) {
        best_x = x;
        best_y = y;
        best_z = z;
        best_side = side;
        best_capacity = capacity;
      }
    }
  }

  std::vector<Vec3> cells;
  cells.reserve(static_cast<std::size_t>(best_capacity));
  const double width = 2.0 / std::sqrt(3.0); // normalized cube edge length
  for (int x = 0; x < best_x; ++x) for (int y = 0; y < best_y; ++y) for (int z = 0; z < best_z; ++z) {
    cells.push_back({
      (static_cast<double>(x) - 0.5 * static_cast<double>(best_x - 1)) * width,
      (static_cast<double>(y) - 0.5 * static_cast<double>(best_y - 1)) * width,
      (static_cast<double>(z) - 0.5 * static_cast<double>(best_z - 1)) * width
    });
  }
  // Spread unused cells throughout the box rather than leaving a conspicuous
  // empty slab. Different workers consequently explore different contact graphs.
  std::shuffle(cells.begin(), cells.end(), rng);

  Packing p;
  p.scale = static_cast<double>(best_side);
  p.poses.resize(static_cast<std::size_t>(cfg.count));
  for (int i = 0; i < cfg.count; ++i) p.poses[static_cast<std::size_t>(i)].position = cells[static_cast<std::size_t>(i)];
  return p;
}

struct LatticePiece {
  Quat rotation;
  std::vector<Vec3> vertices;
  Vec3 width;
  Vec3 box_center;
};

enum class RoundLattice { fcc, bcc };

static LatticePiece lattice_piece(const Polyhedron& piece, const Quat& rotation) {
  LatticePiece result;
  result.rotation = rotation;
  Vec3 lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  Vec3 hi{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  result.vertices.reserve(piece.vertices.size());
  for (const auto& vertex : piece.vertices) {
    const Vec3 v = rotate(rotation, vertex);
    result.vertices.push_back(v);
    lo.x = std::min(lo.x, v.x); lo.y = std::min(lo.y, v.y); lo.z = std::min(lo.z, v.z);
    hi.x = std::max(hi.x, v.x); hi.y = std::max(hi.y, v.y); hi.z = std::max(hi.z, v.z);
  }
  result.width = hi - lo;
  result.box_center = (hi + lo) * 0.5;
  return result;
}

static std::vector<Vec3> contained_lattice_cells(
  const LatticePiece& piece,
  const Polyhedron& shell,
  double scale,
  const Vec3& phase
) {
  const int limit_x = static_cast<int>(std::ceil((scale + 1.0) / piece.width.x)) + 1;
  const int limit_y = static_cast<int>(std::ceil((scale + 1.0) / piece.width.y)) + 1;
  const int limit_z = static_cast<int>(std::ceil((scale + 1.0) / piece.width.z)) + 1;
  std::vector<Vec3> cells;
  for (int x = -limit_x; x <= limit_x; ++x) for (int y = -limit_y; y <= limit_y; ++y) for (int z = -limit_z; z <= limit_z; ++z) {
    const Vec3 box_center{
      (static_cast<double>(x) + phase.x) * piece.width.x,
      (static_cast<double>(y) + phase.y) * piece.width.y,
      (static_cast<double>(z) + phase.z) * piece.width.z
    };
    const Vec3 position = box_center - piece.box_center;
    bool inside = true;
    for (const auto& face : shell.faces) {
      for (const auto& vertex : piece.vertices) {
        if (dot(face.normal, vertex + position) > scale * face.d + 1e-12) {
          inside = false;
          break;
        }
      }
      if (!inside) break;
    }
    if (inside) cells.push_back(position);
  }
  return cells;
}

static bool round_lattice_point(RoundLattice lattice, int x, int y, int z) {
  if (lattice == RoundLattice::fcc) return ((x + y + z) & 1) == 0;
  const int px = std::abs(x) & 1;
  return (std::abs(y) & 1) == px && (std::abs(z) & 1) == px;
}

static double round_lattice_spacing(const LatticePiece& piece, RoundLattice lattice) {
  double spacing = 0.0;
  // Check several neighbour shells. For a convex body, separation by the plane
  // normal to every lattice displacement is conservative; farther neighbours
  // rapidly produce weaker bounds but are included for elongated custom hulls.
  for (int x = -3; x <= 3; ++x) for (int y = -3; y <= 3; ++y) for (int z = -3; z <= 3; ++z) {
    if ((x == 0 && y == 0 && z == 0) || !round_lattice_point(lattice, x, y, z)) continue;
    const Vec3 delta{static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)};
    const Vec3 axis = normalized(delta);
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (const auto& vertex : piece.vertices) {
      const double projection = dot(vertex, axis);
      lo = std::min(lo, projection);
      hi = std::max(hi, projection);
    }
    spacing = std::max(spacing, (hi - lo) / norm(delta));
  }
  return spacing;
}

static std::vector<Vec3> contained_round_lattice_cells(
  const LatticePiece& piece,
  const Polyhedron& shell,
  double scale,
  const Vec3& phase,
  RoundLattice lattice,
  double spacing
) {
  const int limit = static_cast<int>(std::ceil((scale + 1.0) / spacing)) + 2;
  std::vector<Vec3> cells;
  for (int x = -limit; x <= limit; ++x) for (int y = -limit; y <= limit; ++y) for (int z = -limit; z <= limit; ++z) {
    if (!round_lattice_point(lattice, x, y, z)) continue;
    const Vec3 position{
      (static_cast<double>(x) + phase.x) * spacing,
      (static_cast<double>(y) + phase.y) * spacing,
      (static_cast<double>(z) + phase.z) * spacing
    };
    bool inside = true;
    for (const auto& face : shell.faces) {
      for (const auto& vertex : piece.vertices) {
        if (dot(face.normal, vertex + position) > scale * face.d + 1e-12) { inside = false; break; }
      }
      if (!inside) break;
    }
    if (inside) cells.push_back(position);
  }
  return cells;
}

static Packing piece_in_shell_lattice_initial(
  const Polyhedron& piece,
  const Polyhedron& shell,
  const SearchConfig& cfg,
  std::mt19937_64& rng
) {
  double best_scale = std::numeric_limits<double>::infinity();
  std::vector<Vec3> best_cells;
  Quat best_rotation{};

  // A common orientation makes AABB lattice cells a conservative separating
  // construction for every convex piece. Include several orientations because
  // the hull's input axes are often a poor match for sloping shell faces.
  std::vector<Quat> rotations(1);
  for (int i = 0; i < 3; ++i) rotations.push_back(random_q(rng));
  for (const auto& rotation : rotations) {
    const LatticePiece lattice = lattice_piece(piece, rotation);
    for (int px : {0, 1}) for (int py : {0, 1}) for (int pz : {0, 1}) {
      const Vec3 phase{0.5 * static_cast<double>(px), 0.5 * static_cast<double>(py), 0.5 * static_cast<double>(pz)};
      double lo = 0.0;
      double hi = std::max(1.0, 1.35 * std::cbrt(static_cast<double>(cfg.count)));
      while (contained_lattice_cells(lattice, shell, hi, phase).size() < static_cast<std::size_t>(cfg.count)) hi *= 1.5;
      for (int step = 0; step < 40; ++step) {
        const double mid = 0.5 * (lo + hi);
        if (contained_lattice_cells(lattice, shell, mid, phase).size() >= static_cast<std::size_t>(cfg.count)) hi = mid;
        else lo = mid;
      }
      auto cells = contained_lattice_cells(lattice, shell, hi, phase);
      if (hi < best_scale) {
        best_scale = hi;
        best_cells = std::move(cells);
        best_rotation = rotation;
      }
    }

    // Round pieces often waste the corners of their AABB cells. Compete the
    // rectangular seed against face- and body-centred sphere-like lattices.
    for (const RoundLattice pattern : {RoundLattice::fcc, RoundLattice::bcc}) {
      const double spacing = round_lattice_spacing(lattice, pattern);
      for (int px : {0, 1}) for (int py : {0, 1}) for (int pz : {0, 1}) {
        const Vec3 phase{0.5 * static_cast<double>(px), 0.5 * static_cast<double>(py), 0.5 * static_cast<double>(pz)};
        double lo = 0.0;
        double hi = std::max(1.0, 1.35 * std::cbrt(static_cast<double>(cfg.count)));
        while (contained_round_lattice_cells(lattice, shell, hi, phase, pattern, spacing).size() < static_cast<std::size_t>(cfg.count)) hi *= 1.5;
        for (int step = 0; step < 36; ++step) {
          const double mid = 0.5 * (lo + hi);
          if (contained_round_lattice_cells(lattice, shell, mid, phase, pattern, spacing).size() >= static_cast<std::size_t>(cfg.count)) hi = mid;
          else lo = mid;
        }
        auto cells = contained_round_lattice_cells(lattice, shell, hi, phase, pattern, spacing);
        if (hi < best_scale) {
          best_scale = hi;
          best_cells = std::move(cells);
          best_rotation = rotation;
        }
      }
    }
  }

  std::shuffle(best_cells.begin(), best_cells.end(), rng);
  Packing p;
  p.scale = best_scale;
  p.poses.resize(static_cast<std::size_t>(cfg.count));
  for (int i = 0; i < cfg.count; ++i) {
    p.poses[static_cast<std::size_t>(i)].position = best_cells[static_cast<std::size_t>(i)];
    p.poses[static_cast<std::size_t>(i)].rotation = best_rotation;
  }
  return p;
}

static double primary_score(const Packing& p, double penalty) {
  return p.scale + penalty * p.metrics.violation;
}

static double annealing_delta(const Packing& current, const Packing& proposal, double penalty) {
  const double primary = primary_score(proposal, penalty) - primary_score(current, penalty);
  if (std::abs(primary) > 1e-12) return primary;
  // Strictly secondary: hull volume influences only moves tied on shell scale
  // and constraint violation, never trading a larger target shell for compactness.
  const double normalization = std::max(1.0, current.scale * current.scale * current.scale);
  return 0.1 * (proposal.metrics.hull_volume - current.metrics.hull_volume) / normalization;
}

static void order_crystal_core_first(Packing& packing) {
  // Pieces are identical, so reordering poses has no geometric effect. Keeping
  // the outer pieces at the end gives the annealer a stable, seed-independent
  // way to treat a small boundary layer as freely mobile.
  std::sort(packing.poses.begin(), packing.poses.end(), [](const Pose& a, const Pose& b) {
    return norm2(a.position) < norm2(b.position);
  });
}

static Packing one_search(
  const Polyhedron& piece,
  const Polyhedron& shell,
  const SearchConfig& cfg,
  std::uint64_t seed,
  const std::chrono::steady_clock::time_point deadline,
  std::atomic<std::uint64_t>& total_iterations,
  std::function<void(const Packing&)> publish
) {
  std::mt19937_64 rng(seed);
  const double heuristic = std::max(1.25, 1.35 * std::cbrt(static_cast<double>(cfg.count)) * piece.radius / shell.radius);
  const double start_scale = cfg.start_scale > 0.0 ? cfg.start_scale : heuristic;
  const bool lattice_seed = cfg.clearance == 0.0;
  Packing current;
  if (lattice_seed && piece.name == "cube" && shell.name == "cube") current = cube_grid_initial(cfg, rng);
  else if (lattice_seed) current = piece_in_shell_lattice_initial(piece, shell, cfg, rng);
  else current = random_initial(shell, cfg, start_scale, rng);
  if (lattice_seed) order_crystal_core_first(current);
  current.metrics = evaluate(piece, shell, current, cfg.clearance);

  Packing best_any = current;
  Packing best_feasible;
  bool have_feasible = false;

  // Lattice seeds are already in contact, so begin with fine jiggling rather
  // than jumps on the scale of the whole container. The acceptance controller
  // below grows these steps if there is unexpectedly plenty of free space.
  double trans_sigma = 0.045 * piece.radius;
  double rot_sigma = 0.07;
  double scale_sigma = 0.004 * start_scale;
  double temperature = 0.035;
  double penalty = 200.0;
  std::uint64_t iter = 0;
  std::uint64_t window_attempts = 0;
  std::uint64_t window_accepted = 0;
  std::uint64_t last_improvement = 0;
  const int free_piece_count = lattice_seed ? std::min(cfg.count, std::min(6, std::max(2, cfg.count / 5))) : cfg.count;
  const int free_piece_begin = cfg.count - free_piece_count;
  auto random_piece = [&](bool prefer_free) {
    if (prefer_free && rand01(rng) < 0.72) {
      return free_piece_begin + static_cast<int>(rng() % static_cast<std::uint64_t>(free_piece_count));
    }
    return static_cast<int>(rng() % static_cast<std::uint64_t>(cfg.count));
  };
  auto t0 = std::chrono::steady_clock::now();

  // Do not wait for a random move before retaining a constructive seed.  This
  // also guarantees useful output for very short interactive searches.
  if (current.metrics.max_violation <= cfg.tolerance) {
    best_feasible = current;
    have_feasible = true;
    best_feasible.feasible = true;
    publish(best_feasible);
  }

  while (std::chrono::steady_clock::now() < deadline) {
    ++iter;
    ++total_iterations;
    if ((iter % 4000) == 0) {
      // Reheat and tighten constraint pressure. This avoids permanently freezing into one contact graph.
      temperature = std::max(0.015, temperature * 0.85);
      penalty = std::min(2.0e6, penalty * 1.7);
      trans_sigma = std::max(0.008, trans_sigma * 0.86);
      rot_sigma = std::max(0.01, rot_sigma * 0.88);
      scale_sigma = std::max(0.001, scale_sigma * 0.90);
      if (rand01(rng) < 0.22) temperature *= 2.5;
    }

    if ((iter % 500) == 0 && window_attempts != 0) {
      const double acceptance = static_cast<double>(window_accepted) / static_cast<double>(window_attempts);
      const double factor = acceptance < 0.16 ? 0.78 : (acceptance > 0.48 ? 1.18 : 1.0);
      trans_sigma = std::clamp(trans_sigma * factor, 0.0015 * piece.radius, 0.18 * piece.radius);
      rot_sigma = std::clamp(rot_sigma * factor, 0.002, 0.28);
      scale_sigma = std::clamp(scale_sigma * factor, 0.0002 * start_scale, 0.02 * start_scale);
      window_attempts = 0;
      window_accepted = 0;
    }

    // If a contact topology has stopped improving, return to the best known
    // construction and anneal a slightly compressed copy. This is a random
    // restart near useful geometry, not a return to an overlapping cloud.
    if (have_feasible && iter - last_improvement > 6000) {
      current = best_feasible;
      current.scale = std::max(0.25, current.scale * (0.998 + 0.001 * rand01(rng)));
      current.metrics = evaluate(piece, shell, current, cfg.clearance);
      temperature = std::max(temperature, 0.07);
      last_improvement = iter;
    }

    Packing proposal = current;
    const double move = rand01(rng);
    if (move < 0.76) {
      const int i = random_piece(lattice_seed);
      auto& pose = proposal.poses[static_cast<std::size_t>(i)];
      const double mobility = lattice_seed && i >= free_piece_begin ? 2.4 : 1.0;
      if (rand01(rng) < 0.62) {
        pose.position += Vec3{normal01(rng),normal01(rng),normal01(rng)} * (mobility * trans_sigma);
      } else {
        const Quat dq = axis_angle(random_unit(rng), normal01(rng) * (mobility * rot_sigma));
        pose.rotation = normalized(dq * pose.rotation);
      }
    } else if (move < 0.90) {
      // Jiggle a small random neighbourhood together. Multi-piece proposals
      // can cross barriers that reject every intermediate one-piece move.
      const int moved = std::min(cfg.count, 2 + static_cast<int>(rng() % 5));
      for (int k = 0; k < moved; ++k) {
        // Always include a mobile boundary piece, then mix in crystal-core
        // pieces so the interface can reconstruct instead of merely rattling.
        const int i = random_piece(lattice_seed && k == 0);
        auto& pose = proposal.poses[static_cast<std::size_t>(i)];
        const double mobility = lattice_seed && i >= free_piece_begin ? 1.8 : 1.0;
        pose.position += Vec3{normal01(rng),normal01(rng),normal01(rng)} * (0.45 * mobility * trans_sigma);
        if (rand01(rng) < 0.35) {
          const Quat dq = axis_angle(random_unit(rng), normal01(rng) * (0.45 * mobility * rot_sigma));
          pose.rotation = normalized(dq * pose.rotation);
        }
      }
    } else {
      const double shrink_bias = have_feasible ? -0.35*scale_sigma : 0.0;
      proposal.scale = std::max(0.25, proposal.scale + shrink_bias + normal01(rng)*scale_sigma);
    }

    proposal.metrics = evaluate(piece, shell, proposal, cfg.clearance);
    const double delta = annealing_delta(current, proposal, penalty);
    const bool accept = delta <= 0.0 || rand01(rng) < std::exp(-delta/std::max(1e-9,temperature));
    ++window_attempts;
    if (accept) {
      current = std::move(proposal);
      ++window_accepted;
    }

    if (current.metrics.violation < best_any.metrics.violation ||
        (current.metrics.violation <= best_any.metrics.violation*1.000001 && current.scale < best_any.scale)) {
      best_any = current;
    }

    const bool feasible = current.metrics.max_violation <= cfg.tolerance;
    const bool tighter = !have_feasible || current.scale < best_feasible.scale - 1e-12;
    const bool same_shell_more_compact = have_feasible && std::abs(current.scale - best_feasible.scale) <= 1e-12 &&
      current.metrics.hull_volume < best_feasible.metrics.hull_volume;
    if (feasible && (tighter || same_shell_more_compact)) {
      best_feasible = current;
      have_feasible = true;
      best_feasible.feasible = true;
      best_feasible.iterations = iter;
      best_feasible.elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
      publish(best_feasible);
      last_improvement = iter;
      // Once feasible, push gently on the walls immediately.
      current.scale = std::max(0.25, current.scale * 0.9985);
      current.metrics = evaluate(piece, shell, current, cfg.clearance);
    }

    // Rare large shake: move the worst-looking piece to a new region without discarding the whole contact graph.
    if ((iter % 2500) == 0 && rand01(rng) < (lattice_seed ? 0.45 : 0.18)) {
      const int i = random_piece(lattice_seed);
      current.poses[static_cast<std::size_t>(i)].position = random_center_in_shell(shell, current.scale*0.85, rng);
      current.poses[static_cast<std::size_t>(i)].rotation = random_q(rng);
      current.metrics = evaluate(piece, shell, current, cfg.clearance);
    }
  }

  Packing out = have_feasible ? best_feasible : best_any;
  out.feasible = have_feasible;
  out.iterations = iter;
  out.elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  return out;
}

Packing search(const Polyhedron& piece, const Polyhedron& shell, const SearchConfig& config, ProgressCallback progress) {
  if (config.count < 1) throw std::runtime_error("count must be positive");
  if (config.seconds <= 0.0) throw std::runtime_error("seconds must be positive");
  const int thread_count = std::max(1, config.threads);
  const auto deadline = std::chrono::steady_clock::now() +
    std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.seconds));
  std::mutex best_mutex;
  Packing global_best;
  bool have_global = false;
  std::atomic<std::uint64_t> total_iterations{0};

  auto publish = [&](const Packing& candidate) {
    std::lock_guard lock(best_mutex);
    if (!have_global || (candidate.feasible && (!global_best.feasible || candidate.scale < global_best.scale - 1e-12 ||
          (std::abs(candidate.scale - global_best.scale) <= 1e-12 && candidate.metrics.hull_volume < global_best.metrics.hull_volume))) ||
        (!candidate.feasible && !global_best.feasible && candidate.metrics.violation < global_best.metrics.violation)) {
      global_best = candidate;
      have_global = true;
      if (progress) progress(global_best);
    }
  };

  std::vector<std::thread> threads;
  for (int t=0;t<thread_count;t++) {
    threads.emplace_back([&,t] {
      const std::uint64_t seed = config.seed + 0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(t+1);
      Packing local = one_search(piece,shell,config,seed,deadline,total_iterations,publish);
      publish(local);
    });
  }
  for (auto& th : threads) th.join();

  if (!have_global) throw std::runtime_error("search produced no candidate");
  global_best.iterations = total_iterations.load();
  global_best.elapsed = config.seconds;
  global_best.metrics = evaluate(piece,shell,global_best,config.clearance);
  global_best.feasible = global_best.metrics.max_violation <= config.tolerance;
  return global_best;
}

static void write_polyhedron_json(std::ostream& out, const Polyhedron& p) {
  out << "{\n      \"name\": \"" << json_escape(p.name) << "\",\n      \"vertices\": [";
  for (std::size_t i=0;i<p.vertices.size();i++) {
    if (i) out << ',';
    const auto& v=p.vertices[i];
    out << '[' << v.x << ',' << v.y << ',' << v.z << ']';
  }
  out << "],\n      \"faces\": [";
  for (std::size_t i=0;i<p.faces.size();i++) {
    if (i) out << ',';
    out << '[';
    for (std::size_t j=0;j<p.faces[i].vertices.size();j++) { if (j) out << ','; out << p.faces[i].vertices[j]; }
    out << ']';
  }
  out << "]\n    }";
}

void write_result_json(
  const std::string& path,
  const Polyhedron& piece,
  const Polyhedron& shell,
  const SearchConfig& config,
  const Packing& result,
  const std::string& note
) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write result: " + path);
  out << std::setprecision(17);
  out << "{\n  \"format\": \"polyjam-result-v1\",\n";
  out << "  \"piece\": "; write_polyhedron_json(out,piece); out << ",\n";
  out << "  \"shell\": "; write_polyhedron_json(out,shell); out << ",\n";
  out << "  \"count\": " << config.count << ",\n";
  out << "  \"scale\": " << result.scale << ",\n";
  out << "  \"feasible\": " << (result.feasible?"true":"false") << ",\n";
  out << "  \"tolerance\": " << config.tolerance << ",\n";
  out << "  \"clearance\": " << config.clearance << ",\n";
  out << "  \"metrics\": {\"violation\":" << result.metrics.violation
      << ",\"maxViolation\":" << result.metrics.max_violation
      << ",\"containmentSq\":" << result.metrics.containment_sq
      << ",\"overlapSq\":" << result.metrics.overlap_sq
      << ",\"hullVolume\":" << result.metrics.hull_volume
      << ",\"overlappingPairs\":" << result.metrics.overlapping_pairs << "},\n";
  out << "  \"search\": {\"seconds\":" << config.seconds << ",\"threads\":" << config.threads
      << ",\"seed\":" << config.seed << ",\"iterations\":" << result.iterations << "},\n";
  out << "  \"note\": \"" << json_escape(note) << "\",\n";
  out << "  \"poses\": [\n";
  for (std::size_t i=0;i<result.poses.size();i++) {
    const auto& p=result.poses[i];
    out << "    {\"p\":[" << p.position.x << ',' << p.position.y << ',' << p.position.z
        << "],\"q\":[" << p.rotation.w << ',' << p.rotation.x << ',' << p.rotation.y << ',' << p.rotation.z << "]}";
    if (i+1<result.poses.size()) out << ',';
    out << '\n';
  }
  out << "  ]\n}\n";
}

std::string progress_json(const Packing& p, const char* type) {
  std::ostringstream out;
  out << std::setprecision(12);
  out << "{\"type\":\"" << type << "\",\"scale\":" << p.scale
      << ",\"feasible\":" << (p.feasible?"true":"false")
      << ",\"maxViolation\":" << p.metrics.max_violation
      << ",\"violation\":" << p.metrics.violation
      << ",\"hullVolume\":" << p.metrics.hull_volume
      << ",\"iterations\":" << p.iterations
      << ",\"elapsed\":" << p.elapsed << "}";
  return out.str();
}

} // namespace polyjam
