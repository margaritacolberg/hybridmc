// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// event.h contains two types of events that can be found in priority queue:
// collisions between particles and cell crossings

#ifndef HYBRIDMC_EVENT_H
#define HYBRIDMC_EVENT_H

#include <cstdint>
#include <iostream>
#include <queue>
#include <variant>

struct CollisionEvent {
  // time until collision of two hard spheres
  double t;
  // indices of colliding hard spheres
  unsigned int i, j;
  // store copy of collision counters during queue in ni, nj; when an event
  // is processed, compare its counters to ni, nj; discard if mismatch, since
  // this means i, j underwent previous collisions and changed trajectories
  uint64_t ni, nj;
};

struct MinNearestEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os, const MinNearestEvent &ev) {
  return os << "MinNearestEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MaxNearestEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os, const MaxNearestEvent &ev) {
  return os << "MaxNearestEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MinNNearestEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os, const MinNNearestEvent &ev) {
  return os << "MinNNearestEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MaxNNearestEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os, const MaxNNearestEvent &ev) {
  return os << "MaxNNearestEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MinNonlocalInnerEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os,
                                const MinNonlocalInnerEvent &ev) {
  return os << "MinNonlocalInnerEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MaxNonlocalInnerEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os,
                                const MaxNonlocalInnerEvent &ev) {
  return os << "MaxNonlocalInnerEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct MaxNonlocalOuterEvent : public CollisionEvent {};

inline std::ostream &operator<<(std::ostream &os,
                                const MaxNonlocalOuterEvent &ev) {
  return os << "MaxNonlocalOuterEvent{" << ev.t << ", " << ev.i << ", " << ev.j
            << ", " << ev.ni << ", " << ev.nj << "}";
}

struct CellEvent {
  // time until cell crossing of hard sphere
  double t;
  // index of cell crossing hard sphere
  unsigned int i;
  // cell crossing counter of cell crossing hard sphere
  uint64_t ni;
  // 3D index of new cell
  unsigned int ixn, iyn, izn;
  // direction hard sphere is travelling
  enum Wall { xpos, xneg, ypos, yneg, zpos, zneg } wall;
};

struct BeadCellEvent : public CellEvent {};

inline std::ostream &operator<<(std::ostream &os, const BeadCellEvent &ev) {
  return os << "BeadCellEvent{" << ev.t << ", " << ev.i << ", " << ev.ni
            << ", (" << ev.ixn << ", " << ev.iyn << ", " << ev.izn << "), "
            << ev.wall << "}";
}

using Event =
    std::variant<MinNearestEvent, MaxNearestEvent, MinNNearestEvent,
                 MaxNNearestEvent, MinNonlocalInnerEvent, MaxNonlocalInnerEvent,
                 MaxNonlocalOuterEvent, BeadCellEvent>;

// for std::greater<Event>
inline bool operator>(const Event &ev1, const Event &ev2) {
  return std::visit([](auto &&ev1, auto &&ev2) { return ev1.t > ev2.t; }, ev1,
                    ev2);
}

// priority queue is initalized with all possible collisions predicted between
// particles organized from shortest to largest time to collision; as each
// event is processed, queue is updated with new collisions between all
// particle pairs involving a recently-collided particle, and invalid events
// are removed from queue
using EventQueue =
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>;

#endif
