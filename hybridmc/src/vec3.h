// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// vec3.h contains the definition of a 3D vector, and includes commands that
// compare the elements of two 3D vectors to each other

#ifndef HYBRIDMC_VEC3_H
#define HYBRIDMC_VEC3_H

#include <ostream>

struct Vec3 {
  double x, y, z;
  Vec3() = default;
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
};

inline std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
  return os << v.x << " " << v.y << " " << v.z;
}

inline bool operator==(const Vec3 &v1, const Vec3 &v2) {
  return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline bool operator!=(const Vec3 &v1, const Vec3 &v2) {
  return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

#endif
