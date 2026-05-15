from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Vector3d:
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def as_tuple(self):
        return (self.x, self.y, self.z)

    def __add__(self, other):
        return Vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __neg__(self):
        return Vector3d(-self.x, -self.y, -self.z)

    def __mul__(self, value):
        return Vector3d(self.x * value, self.y * value, self.z * value)

    def __truediv__(self, value):
        return Vector3d(self.x / value, self.y / value, self.z / value)

    def norm(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def abs(self):
        return math.sqrt(self.norm())

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        return Vector3d(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def unit(self):
        length = self.abs()
        if length == 0:
            raise ValueError("cannot normalize zero vector")
        return self / length

    def rotated(self, theta, phi, psi):
        cos_t, sin_t = math.cos(theta), math.sin(theta)
        cos_h, sin_h = math.cos(phi), math.sin(phi)
        cos_s, sin_s = math.cos(psi), math.sin(psi)
        ox, oy, oz = self.x, self.y, self.z
        return Vector3d(
            (cos_t * cos_s - sin_t * cos_h * sin_s) * ox
            + (-cos_t * sin_s - sin_t * cos_h * cos_s) * oy
            + (sin_t * sin_h) * oz,
            (sin_t * cos_s + cos_t * cos_h * sin_s) * ox
            + (-sin_t * sin_s + cos_t * cos_h * cos_s) * oy
            + (-cos_t * sin_h) * oz,
            (sin_h * sin_s) * ox + (sin_h * cos_s) * oy + cos_h * oz,
        )

    def axis_rotated(self, axis, theta):
        n = axis.unit()
        b = math.cos(theta)
        c = math.sin(theta)
        a = 1.0 - b
        ox, oy, oz = self.x, self.y, self.z
        return Vector3d(
            (n.x * n.x * a + b) * ox + (n.x * n.y * a - n.z * c) * oy + (n.z * n.x * a + n.y * c) * oz,
            (n.x * n.y * a + n.z * c) * ox + (n.y * n.y * a + b) * oy + (n.y * n.z * a - n.x * c) * oz,
            (n.z * n.x * a - n.y * c) * ox + (n.y * n.z * a + n.x * c) * oy + (n.z * n.z * a + b) * oz,
        )

    def distance_to(self, other):
        return (self - other).abs()


@dataclass
class Point3d:
    x: float
    y: float
    z: float

    def as_tuple(self):
        return (self.x, self.y, self.z)

    def __add__(self, other):
        return Point3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Point3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, value):
        return Point3d(self.x * value, self.y * value, self.z * value)

    def __truediv__(self, value):
        if hasattr(value, "x"):
            return Point3d(self.x / value.x, self.y / value.y, self.z / value.z)
        return Point3d(self.x / value, self.y / value, self.z / value)


def round_half_up(value):
    if value > 0:
        return int(value + 0.5)
    return int(value - 0.5)


def ceili(value):
    return int(math.ceil(value) + 1e-4)
