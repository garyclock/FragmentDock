from __future__ import annotations

from .energy import EnergyCalculator
from .fragment_grid_6d import translation_points_from_config
from .rotation import make_initial_rotations


def dock_full_rotation(conf, receptor, ligands):
    ec = EnergyCalculator(conf.rad_scale)
    rotations = make_initial_rotations(filename=conf.rotangs_file or None)
    translations = sorted(
        translation_points_from_config(conf),
        key=lambda point: (point.x * point.x + point.y * point.y + point.z * point.z, point.x, point.y, point.z),
    )

    scored = []
    max_candidates = max(1, int(conf.poses_per_lig_before_opt))
    for ligand_index, ligand in enumerate(ligands):
        candidate_pairs = []
        for translation in translations:
            for ridx, rotation in enumerate(rotations):
                pose = ligand.rotated_copy(rotation)
                pose.translate(_translation_vec(translation))
                energy = ec.affinity(pose, receptor)
                candidate_pairs.append((pose, energy, ridx, ligand_index))
        candidate_pairs.sort(key=lambda item: (item[1], item[2], item[3]))
        scored.extend(candidate_pairs[:max_candidates])

    scored.sort(key=lambda item: (item[1], item[2], item[3]))
    return [(mol, energy) for mol, energy, _, _ in scored[: max(1, int(conf.poses_per_lig))]]


def _translation_vec(point):
    from .geometry import Vector3d

    return Vector3d(point.x, point.y, point.z)
