from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import LIMIT_ENERGY, XS_TYPE_H
from .fragment_grid import _shifted_grid_values
from .geometry import round_half_up


def _centered_copy(mol):
    centered = mol.copy()
    center = centered.center()
    centered.translate(-center)
    return centered


def _centered_coords(mol):
    heavy = [atom for atom in mol.atoms if atom.xs_type != XS_TYPE_H]
    if not heavy:
        return np.zeros((0, 3), dtype=np.float32)
    coords = np.array([[atom.x, atom.y, atom.z] for atom in heavy], dtype=np.float32)
    coords -= coords.mean(axis=0)
    return coords


def _fragment_geometry_distance(reference, occurrence):
    if reference.size != occurrence.size:
        return float("inf")
    ref = _centered_copy(reference)
    occ = _centered_copy(occurrence)
    total = 0.0
    count = 0
    for atom_a, atom_b in zip(ref.atoms, occ.atoms):
        if atom_a.xs_type == XS_TYPE_H or atom_b.xs_type == XS_TYPE_H:
            continue
        total += (atom_a - atom_b).norm()
        count += 1
    return total / count if count else 0.0


@dataclass(frozen=True)
class FragmentRotationCache:
    coords: np.ndarray


def build_fragment_rotation_cache(normalized_fragment, rotations):
    base = _centered_coords(normalized_fragment)
    rotated = []
    for rotation in rotations:
        frag = normalized_fragment.copy()
        frag.rotate(rotation.x, rotation.y, rotation.z)
        rotated.append(_centered_coords(frag))
    if not rotated:
        return FragmentRotationCache(np.zeros((0,) + base.shape, dtype=np.float32))
    return FragmentRotationCache(np.stack(rotated).astype(np.float32, copy=False))


def _nearest_fragment_rotation_bin_cached(cache, occurrence_fragment):
    if cache.coords.shape[0] == 0:
        return 0
    occurrence = _centered_coords(occurrence_fragment)
    if occurrence.shape != cache.coords.shape[1:]:
        return 0
    diff = cache.coords - occurrence[None, :, :]
    distances = np.sum(diff * diff, axis=(1, 2)) / max(1, occurrence.shape[0])
    return int(np.argmin(distances))


def nearest_fragment_rotation_bin(normalized_fragment, occurrence_fragment, rotations=None):
    if isinstance(normalized_fragment, FragmentRotationCache):
        return _nearest_fragment_rotation_bin_cached(normalized_fragment, occurrence_fragment)
    rotations = list(rotations or [])
    best_idx = 0
    best_dist = float("inf")
    for idx, rotation in enumerate(rotations):
        candidate = normalized_fragment.copy()
        candidate.rotate(rotation.x, rotation.y, rotation.z)
        dist = _fragment_geometry_distance(candidate, occurrence_fragment)
        if dist < best_dist - 1e-12:
            best_idx = idx
            best_dist = dist
    return best_idx


class FragmentInterEnergyGrid6D:
    def __init__(self, orig_frag, rot_angles, atom_grids, distance_grid):
        if not atom_grids:
            raise ValueError("atom_grids is empty")
        self.frag_idx = getattr(orig_frag, "idx", -1)
        self.center = atom_grids[0].center
        self.pitch = atom_grids[0].pitch
        self.num = atom_grids[0].num
        rotations = list(rot_angles)
        if orig_frag.size <= 1:
            rotations = rotations[:1] or []
        self.rotations = rotations
        shape = (int(self.num.x), int(self.num.y), int(self.num.z))
        self.values = np.full((len(rotations),) + shape, LIMIT_ENERGY, dtype=np.float32)
        distance_values = distance_grid.values3d()
        radius = orig_frag.radius()

        for rot_id, rotation in enumerate(rotations):
            frag = orig_frag.copy()
            frag.rotate(rotation.x, rotation.y, rotation.z)
            total = np.zeros(shape, dtype=np.float32)
            for atom in frag.atoms:
                if atom.xs_type == XS_TYPE_H:
                    continue
                dx = round_half_up(atom.x / atom_grids[0].pitch.x)
                dy = round_half_up(atom.y / atom_grids[0].pitch.y)
                dz = round_half_up(atom.z / atom_grids[0].pitch.z)
                total += _shifted_grid_values(atom_grids[atom.xs_type], dx, dy, dz, shape)
                np.minimum(total, LIMIT_ENERGY, out=total)
            if rot_id & 7:
                mask = (distance_values >= 2.0) & (distance_values <= radius + 6.0)
                total = np.where(mask, total, LIMIT_ENERGY).astype(np.float32, copy=False)
            self.values[rot_id] = total

    @property
    def rotation_count(self):
        return self.values.shape[0]

    def get_inter_energy(self, rotation_idx, x, y, z):
        if self.rotation_count == 0:
            return LIMIT_ENERGY
        ridx = min(max(int(rotation_idx), 0), self.rotation_count - 1)
        if x < 0 or y < 0 or z < 0 or x >= self.num.x or y >= self.num.y or z >= self.num.z:
            return LIMIT_ENERGY
        return float(self.values[ridx, int(x), int(y), int(z)])

    def sample_for_search(self, rotation_idx, start, offset, ratio, search_num):
        shape = (int(search_num.x), int(search_num.y), int(search_num.z))
        out = np.full(shape, LIMIT_ENERGY, dtype=np.float32)
        if self.rotation_count == 0:
            return out
        values = self.values[min(max(int(rotation_idx), 0), self.rotation_count - 1)]
        x_idx = int(start.x) + int(offset.x) + np.arange(shape[0]) * int(ratio.x)
        y_idx = int(start.y) + int(offset.y) + np.arange(shape[1]) * int(ratio.y)
        z_idx = int(start.z) + int(offset.z) + np.arange(shape[2]) * int(ratio.z)
        valid_x = (x_idx >= 0) & (x_idx < values.shape[0])
        valid_y = (y_idx >= 0) & (y_idx < values.shape[1])
        valid_z = (z_idx >= 0) & (z_idx < values.shape[2])
        if not (valid_x.any() and valid_y.any() and valid_z.any()):
            return out
        out[np.ix_(valid_x, valid_y, valid_z)] = values[np.ix_(x_idx[valid_x], y_idx[valid_y], z_idx[valid_z])]
        return out
