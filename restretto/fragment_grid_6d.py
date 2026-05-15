from __future__ import annotations

from dataclasses import dataclass

from .geometry import Point3d, Vector3d
from .model import Atom
from .rotation import nearest_rotation_bin


@dataclass(frozen=True)
class Pose6DScore:
    translation_idx: int
    rotation_idx: int
    energy: float


class Fragment6DGrid:
    def __init__(self, translations, rotations, scores):
        self.translations = list(translations)
        self.rotations = list(rotations)
        self.scores = dict(scores)

    def score_pose(self, translation_idx, rotation_idx):
        return self.scores[(translation_idx, rotation_idx)]

    def query_nearest(self, position, rotation):
        translation_idx = min(
            range(len(self.translations)),
            key=lambda idx: (
                _translation_distance(position, self.translations[idx]),
                self.translations[idx].x,
                self.translations[idx].y,
                self.translations[idx].z,
            ),
        )
        energies = [self.score_pose(translation_idx, idx) for idx in range(len(self.rotations))]
        rotation_idx = nearest_rotation_bin(rotation, self.rotations, energies)
        return Pose6DScore(translation_idx, rotation_idx, self.score_pose(translation_idx, rotation_idx))


def _translation_distance(a, b):
    return ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2) ** 0.5


def _translation_vec(point):
    return Vector3d(point.x, point.y, point.z)


def build_fragment_6d_grid(fragment, receptor, translations, rotations, energy_calculator):
    scores = {}
    for tidx, translation in enumerate(translations):
        shift = _translation_vec(translation)
        for ridx, rotation in enumerate(rotations):
            pose = fragment.rotated_copy(rotation)
            pose.translate(shift)
            scores[(tidx, ridx)] = energy_calculator.get_energy(pose, receptor)
    return Fragment6DGrid(translations, rotations, scores)


def translation_points_from_config(conf):
    points = []
    gx = _axis_points(conf.grid.inner_width.x, conf.grid.search_pitch.x)
    gy = _axis_points(conf.grid.inner_width.y, conf.grid.search_pitch.y)
    gz = _axis_points(conf.grid.inner_width.z, conf.grid.search_pitch.z)
    for x in gx:
        for y in gy:
            for z in gz:
                points.append(Point3d(x, y, z))
    return points


def _axis_points(width, pitch):
    half = width / 2.0
    count = int(round(width / pitch))
    return [-half + idx * pitch for idx in range(count + 1)]
