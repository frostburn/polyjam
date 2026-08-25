#include "geometry.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace polyjam {

Vec3& Vec3::operator+=(const Vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
Vec3& Vec3::operator-=(const Vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
Vec3& Vec3::operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
Vec3& Vec3::operator/=(double s) { x /= s; y /= s; z /= s; return *this; }
Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
Vec3 operator*(Vec3 a, double s) { return a *= s; }
Vec3 operator*(double s, Vec3 a) { return a *= s; }
Vec3 operator/(Vec3 a, double s) { return a /= s; }
double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
double norm2(const Vec3& v) { return dot(v, v); }
double norm(const Vec3& v) { return std::sqrt(norm2(v)); }
Vec3 normalized(const Vec3& v) {
  const double n = norm(v);
  if (n == 0.0) return {};
  return v / n;
}

Quat normalized(const Quat& q) {
  const double n = std::sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
  if (n == 0.0) return {};
  return {q.w/n, q.x/n, q.y/n, q.z/n};
}
Quat operator*(const Quat& a, const Quat& b) {
  return {
    a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
    a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
    a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
    a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
  };
}
Quat axis_angle(const Vec3& axis, double angle) {
  const Vec3 u = normalized(axis);
  const double h = 0.5 * angle;
  const double s = std::sin(h);
  return normalized({std::cos(h), u.x*s, u.y*s, u.z*s});
}
Vec3 rotate(const Quat& q0, const Vec3& v) {
  const Quat q = normalized(q0);
  const Vec3 u{q.x, q.y, q.z};
  const double s = q.w;
  return 2.0 * dot(u, v) * u + (s*s - dot(u,u)) * v + 2.0*s*cross(u, v);
}

static std::uint64_t splitmix64(std::uint64_t& x) {
  std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
static double uniform01(std::uint64_t& state) {
  return static_cast<double>(splitmix64(state) >> 11) * (1.0 / 9007199254740992.0);
}
Quat random_quat(std::uint64_t& state) {
  // Shoemake: uniform orientation on S^3.
  const double u1 = uniform01(state);
  const double u2 = uniform01(state);
  const double u3 = uniform01(state);
  const double a = std::sqrt(1.0-u1);
  const double b = std::sqrt(u1);
  const double t1 = 2.0 * std::numbers::pi * u2;
  const double t2 = 2.0 * std::numbers::pi * u3;
  return normalized({b*std::cos(t2), a*std::sin(t1), a*std::cos(t1), b*std::sin(t2)});
}

static std::vector<Vec3> cube_vertices() {
  std::vector<Vec3> v;
  for (int x : {-1, 1}) for (int y : {-1, 1}) for (int z : {-1, 1}) v.push_back({double(x),double(y),double(z)});
  return v;
}
static std::vector<Vec3> tetra_vertices() {
  return {{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
}
static std::vector<Vec3> octa_vertices() {
  return {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
}
static std::vector<Vec3> icosa_vertices() {
  const double p = (1.0 + std::sqrt(5.0)) / 2.0;
  std::vector<Vec3> v;
  for (int a : {-1,1}) for (int b : {-1,1}) {
    v.push_back({0.0, double(a), double(b)*p});
    v.push_back({double(a), double(b)*p, 0.0});
    v.push_back({double(b)*p, 0.0, double(a)});
  }
  return v;
}
static std::vector<Vec3> dodeca_vertices() {
  const double p = (1.0 + std::sqrt(5.0)) / 2.0;
  const double q = 1.0 / p;
  std::vector<Vec3> v;
  for (int x : {-1,1}) for (int y : {-1,1}) for (int z : {-1,1}) v.push_back({double(x),double(y),double(z)});
  for (int a : {-1,1}) for (int b : {-1,1}) {
    v.push_back({0.0, double(a)*q, double(b)*p});
    v.push_back({double(a)*q, double(b)*p, 0.0});
    v.push_back({double(b)*p, 0.0, double(a)*q});
  }
  return v;
}

std::vector<std::string> builtin_names() { return {"tetra", "cube", "octa", "icosa", "dodeca"}; }

static void recenter_and_normalize(std::vector<Vec3>& vertices) {
  if (vertices.empty()) throw std::runtime_error("polyhedron has no vertices");
  Vec3 c{};
  for (const auto& v : vertices) c += v;
  c /= static_cast<double>(vertices.size());
  double r = 0.0;
  for (auto& v : vertices) { v -= c; r = std::max(r, norm(v)); }
  if (r <= 1e-12) throw std::runtime_error("degenerate polyhedron");
  for (auto& v : vertices) v /= r;
}

static bool same_plane(const Face& a, const Vec3& n, double d) {
  return dot(a.normal, n) > 1.0 - 1e-7 && std::abs(a.d - d) < 1e-7;
}

Polyhedron make_convex_hull(std::string name, std::vector<Vec3> vertices) {
  recenter_and_normalize(vertices);
  if (vertices.size() < 4) throw std::runtime_error("need at least four vertices for a 3D convex polyhedron");

  std::vector<Face> planes;
  const double eps = 1e-8;
  const int n = static_cast<int>(vertices.size());
  for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) for (int k=j+1;k<n;k++) {
    Vec3 nn = cross(vertices[j]-vertices[i], vertices[k]-vertices[i]);
    const double len = norm(nn);
    if (len < eps) continue;
    nn /= len;
    double d = dot(nn, vertices[i]);
    double max_pos = -std::numeric_limits<double>::infinity();
    double min_pos = std::numeric_limits<double>::infinity();
    for (const auto& v : vertices) {
      const double s = dot(nn, v) - d;
      max_pos = std::max(max_pos, s);
      min_pos = std::min(min_pos, s);
    }
    if (max_pos > eps && min_pos < -eps) continue; // not a support plane
    if (max_pos > eps) { nn = -nn; d = -d; } // make all vertices <= d
    if (d < 0.0) { nn = -nn; d = -d; } // origin should be inside; stable outward orientation
    bool duplicate = false;
    for (const auto& f : planes) if (same_plane(f, nn, d)) { duplicate = true; break; }
    if (!duplicate) planes.push_back({{}, nn, d});
  }
  if (planes.size() < 4) throw std::runtime_error("input vertices do not form a full-dimensional convex hull");

  for (auto& face : planes) {
    for (int i=0;i<n;i++) if (std::abs(dot(face.normal, vertices[i]) - face.d) < 2e-7) face.vertices.push_back(i);
    if (face.vertices.size() < 3) throw std::runtime_error("internal hull error: support plane with <3 vertices");

    Vec3 center{};
    for (int idx : face.vertices) center += vertices[idx];
    center /= static_cast<double>(face.vertices.size());
    Vec3 u = normalized(vertices[face.vertices[0]] - center);
    Vec3 v = normalized(cross(face.normal, u));
    std::sort(face.vertices.begin(), face.vertices.end(), [&](int a, int b) {
      const Vec3 da = vertices[a]-center;
      const Vec3 db = vertices[b]-center;
      const double aa = std::atan2(dot(da,v), dot(da,u));
      const double ab = std::atan2(dot(db,v), dot(db,u));
      return aa < ab;
    });
    // Ensure winding is outward when viewed from outside.
    if (face.vertices.size() >= 3) {
      const Vec3 e1 = vertices[face.vertices[1]] - vertices[face.vertices[0]];
      const Vec3 e2 = vertices[face.vertices[2]] - vertices[face.vertices[1]];
      if (dot(cross(e1,e2), face.normal) < 0.0) std::reverse(face.vertices.begin(), face.vertices.end());
    }
  }

  std::set<std::pair<int,int>> edge_set;
  for (const auto& f : planes) {
    for (std::size_t i=0;i<f.vertices.size();i++) {
      int a = f.vertices[i];
      int b = f.vertices[(i+1)%f.vertices.size()];
      if (a>b) std::swap(a,b);
      edge_set.emplace(a,b);
    }
  }

  Polyhedron p;
  p.name = std::move(name);
  p.vertices = std::move(vertices);
  p.faces = std::move(planes);
  p.edges.assign(edge_set.begin(), edge_set.end());
  p.radius = 1.0;
  return p;
}

Polyhedron builtin_polyhedron(const std::string& name) {
  if (name == "cube") return make_convex_hull(name, cube_vertices());
  if (name == "tetra" || name == "tetrahedron") return make_convex_hull("tetra", tetra_vertices());
  if (name == "octa" || name == "octahedron") return make_convex_hull("octa", octa_vertices());
  if (name == "icosa" || name == "icosahedron") return make_convex_hull("icosa", icosa_vertices());
  if (name == "dodeca" || name == "dodecahedron") return make_convex_hull("dodeca", dodeca_vertices());
  throw std::runtime_error("unknown builtin polyhedron: " + name);
}

Polyhedron load_obj_convex_hull(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open OBJ: " + path);
  std::vector<Vec3> vertices;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string tag;
    ss >> tag;
    if (tag == "v") {
      Vec3 v;
      if (ss >> v.x >> v.y >> v.z) vertices.push_back(v);
    }
  }
  if (vertices.empty()) throw std::runtime_error("OBJ has no 'v' vertices: " + path);
  return make_convex_hull(path, std::move(vertices));
}

std::vector<Vec3> transformed_vertices(const Polyhedron& p, const Pose& pose) {
  std::vector<Vec3> out;
  out.reserve(p.vertices.size());
  for (const auto& v : p.vertices) out.push_back(rotate(pose.rotation, v) + pose.position);
  return out;
}

static std::pair<double,double> project(const std::vector<Vec3>& vertices, const Vec3& axis) {
  double lo = std::numeric_limits<double>::infinity();
  double hi = -std::numeric_limits<double>::infinity();
  for (const auto& v : vertices) {
    const double p = dot(v, axis);
    lo = std::min(lo,p); hi = std::max(hi,p);
  }
  return {lo,hi};
}

SeparationResult sat_test(
  const Polyhedron& shape,
  const Pose& a,
  const std::vector<Vec3>& av,
  const Pose& b,
  const std::vector<Vec3>& bv,
  double axis_epsilon
) {
  double min_overlap = std::numeric_limits<double>::infinity();
  auto test_axis = [&](Vec3 axis) -> bool {
    const double n2 = norm2(axis);
    if (n2 < axis_epsilon*axis_epsilon) return false;
    axis /= std::sqrt(n2);
    const auto [alo,ahi] = project(av, axis);
    const auto [blo,bhi] = project(bv, axis);
    const double overlap = std::min(ahi,bhi) - std::max(alo,blo);
    if (overlap <= 0.0) return true;
    min_overlap = std::min(min_overlap, overlap);
    return false;
  };

  for (const auto& f : shape.faces) if (test_axis(rotate(a.rotation, f.normal))) return {true,0.0};
  for (const auto& f : shape.faces) if (test_axis(rotate(b.rotation, f.normal))) return {true,0.0};
  for (const auto& ea : shape.edges) {
    const Vec3 ad = rotate(a.rotation, shape.vertices[ea.second]-shape.vertices[ea.first]);
    for (const auto& eb : shape.edges) {
      const Vec3 bd = rotate(b.rotation, shape.vertices[eb.second]-shape.vertices[eb.first]);
      if (test_axis(cross(ad,bd))) return {true,0.0};
    }
  }
  if (!std::isfinite(min_overlap)) min_overlap = 0.0;
  return {false,min_overlap};
}

ContainmentResult containment_violation(
  const Polyhedron& shell,
  double shell_scale,
  const std::vector<Vec3>& world_vertices
) {
  ContainmentResult r;
  for (const auto& f : shell.faces) {
    const double limit = shell_scale * f.d;
    for (const auto& v : world_vertices) {
      const double excess = dot(f.normal, v) - limit;
      if (excess > 0.0) {
        r.sum_sq += excess*excess;
        r.max_excess = std::max(r.max_excess, excess);
      }
    }
  }
  return r;
}

std::string json_escape(const std::string& s) {
  std::string out;
  out.reserve(s.size()+8);
  for (char c : s) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out += c; break;
    }
  }
  return out;
}

} // namespace polyjam
