#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace polyjam {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  Vec3& operator+=(const Vec3& rhs);
  Vec3& operator-=(const Vec3& rhs);
  Vec3& operator*=(double s);
  Vec3& operator/=(double s);
};

Vec3 operator+(Vec3 a, const Vec3& b);
Vec3 operator-(Vec3 a, const Vec3& b);
Vec3 operator-(const Vec3& a);
Vec3 operator*(Vec3 a, double s);
Vec3 operator*(double s, Vec3 a);
Vec3 operator/(Vec3 a, double s);

double dot(const Vec3& a, const Vec3& b);
Vec3 cross(const Vec3& a, const Vec3& b);
double norm2(const Vec3& v);
double norm(const Vec3& v);
Vec3 normalized(const Vec3& v);

struct Quat {
  double w = 1.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Quat normalized(const Quat& q);
Quat operator*(const Quat& a, const Quat& b);
Quat axis_angle(const Vec3& axis, double angle);
Vec3 rotate(const Quat& q, const Vec3& v);
Quat random_quat(std::uint64_t& state);

struct Face {
  std::vector<int> vertices;
  Vec3 normal;
  double d = 0.0; // outward plane: normal . x <= d
};

struct Polyhedron {
  std::string name;
  std::vector<Vec3> vertices;
  std::vector<Face> faces;
  std::vector<std::pair<int, int>> edges;
  double radius = 1.0;
};

Polyhedron builtin_polyhedron(const std::string& name);
Polyhedron load_obj_convex_hull(const std::string& path);
Polyhedron make_convex_hull(std::string name, std::vector<Vec3> vertices);
std::vector<std::string> builtin_names();

struct Pose {
  Vec3 position;
  Quat rotation;
};

std::vector<Vec3> transformed_vertices(const Polyhedron& p, const Pose& pose);

struct SeparationResult {
  bool separated = false;
  double penetration = 0.0; // 0 if separated, otherwise min SAT overlap
};

SeparationResult sat_test(
  const Polyhedron& shape,
  const Pose& a,
  const std::vector<Vec3>& av,
  const Pose& b,
  const std::vector<Vec3>& bv,
  double axis_epsilon = 1e-10
);

struct ContainmentResult {
  double sum_sq = 0.0;
  double max_excess = 0.0;
};

ContainmentResult containment_violation(
  const Polyhedron& shell,
  double shell_scale,
  const std::vector<Vec3>& world_vertices
);

std::string json_escape(const std::string& s);

} // namespace polyjam
