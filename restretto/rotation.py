from __future__ import annotations

import math
from pathlib import Path

from .geometry import Vector3d


def _angle(a, b):
    denom = a.abs() * b.abs()
    if denom == 0:
        return 0.0
    value = max(-1.0, min(1.0, a.dot(b) / denom))
    return math.acos(value)


def _get_rot(v1, v2):
    vec1 = Vector3d(v1.x, v1.y, v1.z)
    vec2 = v1.cross(v2)
    psi = _angle(Vector3d(0.0, 1.0, 0.0), Vector3d(vec2.x, vec2.y, 0.0))
    phi = _angle(Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, math.sqrt(vec2.x * vec2.x + vec2.y * vec2.y), vec2.z))
    if vec2.x < 0:
        psi = -psi
    vec1 = vec1.rotated(0.0, phi, psi)
    theta = -_angle(Vector3d(1.0, 0.0, 0.0), vec1)
    if vec1.y < 0:
        theta = -theta
    return Vector3d(-psi, -phi, -theta)


def make_rotations_60():
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    pi = math.acos(-1.0)
    pole_a = Vector3d(0.0, golden * golden * golden, golden * golden).unit()
    pole_b = Vector3d(0.0, 1.0, golden * golden).unit()
    p = _angle(pole_a, pole_b)
    pole_a = Vector3d(1.0, 0.0, 0.0)
    pole_b = Vector3d(math.cos(p), math.sin(p), 0.0)
    order = [0, 0, 1, 1, 1, 2, 1, 2, 2, 2, 3, -1]

    rotations = []
    for o in order:
        for _ in range(5):
            rotations.append(_get_rot(pole_a, pole_b))
            pole_b = pole_b.axis_rotated(pole_a, 2.0 * pi / 5.0)
        pole_b = pole_b.axis_rotated(pole_a, 2.0 * o * pi / 5.0)
        pole_a = pole_a.axis_rotated(pole_b, 4.0 * pi / 3.0)
    return rotations


def _wrap_angle(delta):
    return (delta + math.pi) % (2.0 * math.pi) - math.pi


def rotation_distance(a, b):
    return math.sqrt(
        _wrap_angle(a.x - b.x) ** 2
        + _wrap_angle(a.y - b.y) ** 2
        + _wrap_angle(a.z - b.z) ** 2
    )


def make_initial_rotations(reference_mode=True, filename=None):
    if filename:
        return read_rotations(filename)
    if not reference_mode:
        return [Vector3d(0.0, 0.0, 0.0)]
    return make_rotations_60()


def read_rotations(filename):
    rotations = []
    for raw in Path(filename).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [float(part.strip()) for part in line.split(",")]
        if len(parts) != 3:
            raise ValueError("rotation line must have three comma-separated angles")
        rotations.append(Vector3d(parts[0], parts[1], parts[2]))
    if not rotations:
        raise ValueError("rotation file did not contain any rotations")
    return rotations


def choose_better_rotation(candidate_a, candidate_b):
    idx_a, dist_a, energy_a = candidate_a
    idx_b, dist_b, energy_b = candidate_b
    if dist_b < dist_a - 1e-12:
        return candidate_b
    if dist_a < dist_b - 1e-12:
        return candidate_a
    if energy_b < energy_a - 1e-12:
        return candidate_b
    if energy_a < energy_b - 1e-12:
        return candidate_a
    return candidate_b if idx_b < idx_a else candidate_a


def nearest_rotation_bin(rotation, bins, energies=None):
    if not bins:
        raise ValueError("rotation bins must not be empty")
    energies = energies if energies is not None else [0.0] * len(bins)
    best = (0, rotation_distance(rotation, bins[0]), energies[0])
    for idx, candidate in enumerate(bins[1:], start=1):
        best = choose_better_rotation(best, (idx, rotation_distance(rotation, candidate), energies[idx]))
    return best[0]
