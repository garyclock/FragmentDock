from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

from .config import parse_in_file
from .config import ReuseStrategy
from .constants import (
    LIMIT_ENERGY,
    XS_STRINGS,
    XS_TYPE_DUMMY,
    XS_TYPE_H,
    XS_TYPE_SIZE,
    xs_hbond,
    xs_is_hydrophobic,
    xs_radius,
)
from .decompose import decompose_molecule
from .docking_full_rotation import dock_full_rotation
from .energy import EnergyCalculator, PRECI, SZ, TERM_WEIGHTS, THRESHOLD
from .fragment_grid import (
    FragmentInterEnergyGrid,
    FragmentInterEnergyGridContainer,
    build_offline_schedule,
    make_distance_grid,
)
from .geometry import Point3d, Vector3d, ceili, round_half_up
from .grid import InterEnergyGrid
from .io import _canonical_smiles, _to_obmol, canonical_labels, read_molecules, write_sdf_like
from .optimizer import OptimizerGrid
from .rotation import make_initial_rotations, make_rotations_60


def _resolve_config_path(config_path, value):
    path = Path(value)
    if path.is_absolute():
        return path
    candidate = Path(config_path).parent / path
    if candidate.exists():
        return candidate
    return path


def _load_inputs(config_path):
    conf = parse_in_file(config_path)
    return _load_inputs_from_conf(config_path, conf)


def _load_inputs_from_conf(config_path, conf):
    receptor_path = _resolve_config_path(config_path, conf.receptor_file)
    ligand_paths = [_resolve_config_path(config_path, ligand) for ligand in conf.ligand_files]
    receptor = read_molecules(receptor_path)[0]
    ligands = []
    for path in ligand_paths:
        if path.exists():
            ligands.extend(read_molecules(path))
    return conf, receptor, ligands


def _apply_common_overrides(conf, args):
    if getattr(args, "ligand", None):
        conf.ligand_files = list(args.ligand)
    if getattr(args, "receptor", None):
        conf.receptor_file = args.receptor
    if getattr(args, "output", None):
        conf.output_file = args.output
    if getattr(args, "grid", None):
        conf.grid_folder = args.grid
    if getattr(args, "memsize", None) is not None:
        conf.mem_size = args.memsize
    if getattr(args, "score_only", False):
        conf.score_only = True
    if getattr(args, "local_only", False):
        conf.local_only = True
    if getattr(args, "local_max_rmsd", None) is not None:
        conf.local_max_rmsd = args.local_max_rmsd
    if getattr(args, "log", None):
        conf.log_file = args.log
    if getattr(args, "poses_per_lig", None) is not None:
        conf.poses_per_lig = args.poses_per_lig
    if getattr(args, "min_rmsd", None) is not None:
        conf.pose_min_rmsd = args.min_rmsd
    if getattr(args, "no_local_opt", False):
        conf.no_local_opt = True
    if getattr(args, "poses_per_lig_before_opt", None) is not None:
        conf.poses_per_lig_before_opt = args.poses_per_lig_before_opt
    if getattr(args, "output_score_threshold", None) is not None:
        conf.output_score_threshold = args.output_score_threshold
    if getattr(args, "dxgrid", None):
        conf.dxgrid_folder = args.dxgrid
    conf.check_validity()
    return conf


def atomgrid_gen(args):
    conf = _apply_common_overrides(parse_in_file(args.config), args)
    conf, receptor, _ = _load_inputs_from_conf(args.config, conf)
    grid_dir = Path(conf.grid_folder)
    if not grid_dir.is_absolute():
        grid_dir = Path(args.config).parent / grid_dir
    grid_dir.mkdir(parents=True, exist_ok=True)
    num = Point3d(
        ceili(conf.grid.outer_width.x / 2.0 / conf.grid.score_pitch.x) * 2 + 1,
        ceili(conf.grid.outer_width.y / 2.0 / conf.grid.score_pitch.y) * 2 + 1,
        ceili(conf.grid.outer_width.z / 2.0 / conf.grid.score_pitch.z) * 2 + 1,
    )
    ec = EnergyCalculator(conf.rad_scale)
    for xs_type in range(XS_TYPE_SIZE):
        grid = _build_atom_grid(conf.grid.center, conf.grid.score_pitch, num, xs_type, receptor, ec)
        grid.write_grid(grid_dir / ("%s.grid" % XS_STRINGS[xs_type]))
    return 0


def _energy_lookup_table(t1, t2, rad_scale):
    dist = np.arange(SZ, dtype=np.float64) / float(PRECI)
    d = dist - (xs_radius(t1) + xs_radius(t2)) * rad_scale
    value = TERM_WEIGHTS[2] * np.where(d > 0.0, 0.0, d * d)
    if t1 not in {XS_TYPE_H, XS_TYPE_DUMMY} and t2 not in {XS_TYPE_H, XS_TYPE_DUMMY}:
        value += TERM_WEIGHTS[0] * np.exp(-((d * 2.0) ** 2))
        value += TERM_WEIGHTS[1] * np.exp(-(((d - 3.0) * 0.5) ** 2))
        if xs_is_hydrophobic(t1) and xs_is_hydrophobic(t2):
            value += TERM_WEIGHTS[3] * np.where(d >= 1.5, 0.0, np.where(d <= 0.5, 1.0, 1.5 - d))
        if xs_hbond(t1, t2):
            value += TERM_WEIGHTS[4] * np.where(d >= 0.0, 0.0, np.where(d <= -0.7, 1.0, d * -1.428571))
    return value.astype(np.float32)


def _grid_axis(center, pitch, count):
    return (np.arange(int(count), dtype=np.float64) - (int(count) - 1) / 2.0) * pitch + center


def _build_atom_grid(center, pitch, num, xs_type, receptor, ec):
    grid = InterEnergyGrid(center, pitch, num, 0.0)
    values = grid.values3d()
    xs = _grid_axis(center.x, pitch.x, num.x)
    ys = _grid_axis(center.y, pitch.y, num.y)
    zs = _grid_axis(center.z, pitch.z, num.z)
    radius_sq = float(THRESHOLD * THRESHOLD)
    tables = {}

    for atom in receptor.atoms:
        if atom.xs_type == XS_TYPE_H:
            continue
        center_x = (atom.x - center.x) / pitch.x + (int(num.x) - 1) / 2.0
        center_y = (atom.y - center.y) / pitch.y + (int(num.y) - 1) / 2.0
        center_z = (atom.z - center.z) / pitch.z + (int(num.z) - 1) / 2.0
        rx = THRESHOLD / pitch.x
        ry = THRESHOLD / pitch.y
        rz = THRESHOLD / pitch.z
        x0 = max(0, int(math.floor(center_x - rx)))
        y0 = max(0, int(math.floor(center_y - ry)))
        z0 = max(0, int(math.floor(center_z - rz)))
        x1 = min(int(num.x), int(math.ceil(center_x + rx)) + 1)
        y1 = min(int(num.y), int(math.ceil(center_y + ry)) + 1)
        z1 = min(int(num.z), int(math.ceil(center_z + rz)) + 1)
        if x0 >= x1 or y0 >= y1 or z0 >= z1:
            continue

        table = tables.setdefault(atom.xs_type, _energy_lookup_table(xs_type, atom.xs_type, ec.rad_scale))
        dx2 = (xs[x0:x1] - atom.x) ** 2
        dy2 = (ys[y0:y1] - atom.y) ** 2
        dz2 = (zs[z0:z1] - atom.z) ** 2
        dist_sq = dx2[:, None, None] + dy2[None, :, None] + dz2[None, None, :]
        mask = dist_sq <= radius_sq
        if not np.any(mask):
            continue
        idx = (np.sqrt(dist_sq, dtype=np.float64) * PRECI).astype(np.int32)
        mask &= idx < SZ
        if not np.any(mask):
            continue
        target = values[x0:x1, y0:y1, z0:z1]
        active = target < LIMIT_ENERGY
        mask &= active
        if not np.any(mask):
            continue
        target[mask] += table[idx[mask]]
        target[target >= LIMIT_ENERGY] = LIMIT_ENERGY
    return grid


def score_only(args):
    conf = _apply_common_overrides(parse_in_file(args.config), args)
    _, receptor, ligands = _load_inputs_from_conf(args.config, conf)
    ec = EnergyCalculator()
    for ligand in ligands:
        print("Title: %s" % ligand.title)
        print("Affinity: %.5f" % ec.affinity(ligand, receptor))
    return 0


def intraenergy_only(args):
    conf = _apply_common_overrides(parse_in_file(args.config), args)
    _, _, ligands = _load_inputs_from_conf(args.config, conf)
    ec = EnergyCalculator()
    for ligand in ligands:
        print("Title: %s" % ligand.title)
        print("IntraEnergy: %.5f" % ec.calc_intra_energy(ligand))
    return 0


def decompose(args):
    if getattr(args, "config", None):
        _, _, ligands = _load_inputs(args.config)
    else:
        ligands = []
        for ligand_path in args.ligand or []:
            ligands.extend(read_molecules(ligand_path))
    scored = []
    annotated = []
    for ligand in ligands:
        fragments = decompose_molecule(ligand)
        frag_names = []
        for fragment in fragments:
            frag_names.append(fragment.smiles or fragment.title)
            scored.append((fragment, 0.0))
        annotated_ligand = ligand.copy()
        annotated_ligand.properties = {"fragment_info": ",".join(frag_names)}
        annotated.append((annotated_ligand, 0.0))
    if getattr(args, "fragment", None):
        write_sdf_like(args.fragment, scored)
        write_sdf_like(args.output, annotated)
    else:
        write_sdf_like(args.output, scored)
    return 0


def conformer_docking(args):
    conf = _apply_common_overrides(parse_in_file(args.config), args)
    conf, receptor, ligands = _load_inputs_from_conf(args.config, conf)
    if getattr(args, "full_rotation", False):
        scored = dock_full_rotation(conf, receptor, ligands)
    elif conf.score_only or conf.local_only:
        scored = _dock_with_atom_grids(args.config, conf, ligands)
    else:
        scored = _dock_with_fragment_grids(args.config, conf, receptor, ligands)
    output = Path(conf.output_file)
    if not output.is_absolute():
        output = Path(args.config).parent / output
    write_sdf_like(output, scored[: max(1, int(conf.poses_per_lig))])
    if not (conf.score_only or conf.local_only):
        _write_score_csv(output, scored)
    return 0


def _date_string():
    now = time.localtime()
    return "%02d_%02d_%02d_%02d_%02d" % (now.tm_mon - 1, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)


def _write_score_csv(output, scored):
    csv_path = Path(str(output) + "fraggrid__" + _date_string() + ".csv")
    best_by_identifier = {}
    for mol, score in scored:
        identifier = mol.identifier
        if identifier not in best_by_identifier or score < best_by_identifier[identifier]:
            best_by_identifier[identifier] = score
    with csv_path.open("w", encoding="utf-8") as stream:
        for identifier, score in best_by_identifier.items():
            stream.write("%s,%s\n" % (identifier, score))


def _read_atom_grids(config_path, conf):
    grid_dir = Path(conf.grid_folder)
    if not grid_dir.is_absolute():
        grid_dir = Path(config_path).parent / grid_dir
    grids = []
    for name in XS_STRINGS:
        grids.append(InterEnergyGrid.parse_grid(grid_dir / ("%s.grid" % name)))
    return grids


def _score_on_atom_grids(molecule, atom_grids):
    total = 0.0
    for atom in molecule.atoms:
        if atom.xs_type == XS_TYPE_H:
            continue
        total += atom_grids[atom.xs_type].get_inter_energy(atom)
        if total >= LIMIT_ENERGY:
            return LIMIT_ENERGY
    return total


def _search_grid_num(conf):
    return Point3d(
        ceili(conf.grid.inner_width.x / 2.0 / conf.grid.search_pitch.x) * 2 + 1,
        ceili(conf.grid.inner_width.y / 2.0 / conf.grid.search_pitch.y) * 2 + 1,
        ceili(conf.grid.inner_width.z / 2.0 / conf.grid.search_pitch.z) * 2 + 1,
    )


def _score_grid_num(conf):
    return Point3d(
        ceili(conf.grid.outer_width.x / 2.0 / conf.grid.score_pitch.x) * 2 + 1,
        ceili(conf.grid.outer_width.y / 2.0 / conf.grid.score_pitch.y) * 2 + 1,
        ceili(conf.grid.outer_width.z / 2.0 / conf.grid.score_pitch.z) * 2 + 1,
    )


def _grid_ratio(conf):
    return Point3d(
        round_half_up(conf.grid.search_pitch.x / conf.grid.score_pitch.x),
        round_half_up(conf.grid.search_pitch.y / conf.grid.score_pitch.y),
        round_half_up(conf.grid.search_pitch.z / conf.grid.score_pitch.z),
    )


def _to_score_num(k, score_num, search_num, ratio):
    return Point3d(
        int(score_num.x / 2) + (-int(search_num.x / 2) + k) * int(ratio.x),
        int(score_num.y / 2) + (-int(search_num.y / 2) + k) * int(ratio.y),
        int(score_num.z / 2) + (-int(search_num.z / 2) + k) * int(ratio.z),
    )


def _round_fragment_offset(vec, pitch):
    return Point3d(
        round_half_up(vec.x / pitch.x),
        round_half_up(vec.y / pitch.y),
        round_half_up(vec.z / pitch.z),
    )


def _fragment_signature(fragment):
    atoms = tuple(atom.xs_type for atom in fragment.atoms)
    bonds = tuple(sorted((min(bond.atom_id1, bond.atom_id2), max(bond.atom_id1, bond.atom_id2), bond.is_rotor) for bond in fragment.bonds))
    return atoms, bonds


def _angle(a, b):
    denom = a.abs() * b.abs()
    if denom == 0:
        return 0.0
    value = max(-1.0, min(1.0, a.dot(b) / denom))
    return math.acos(value)


def _fragment_reference_triplet(fragment):
    if fragment.size == 1:
        return 0, -1, -1
    tri = [0, 1, -1]
    atom_by_id = {atom.id: atom for atom in fragment.atoms}
    pi = math.acos(-1.0)
    for atom in fragment.atoms:
        if atom.id >= len(fragment.bond_ids):
            continue
        bond_ids = fragment.bond_ids[atom.id]
        for j in range(len(bond_ids)):
            for k in range(j + 1, len(bond_ids)):
                bond1 = fragment.bonds[bond_ids[j]]
                bond2 = fragment.bonds[bond_ids[k]]
                atom_id1 = bond1.atom_id1 + bond1.atom_id2 - atom.id
                atom_id2 = bond2.atom_id1 + bond2.atom_id2 - atom.id
                if atom_id1 not in atom_by_id or atom_id2 not in atom_by_id:
                    continue
                vec1 = atom_by_id[atom_id1] - atom
                vec2 = atom_by_id[atom_id2] - atom
                angle = _angle(vec1, vec2)
                if angle > pi * 0.1 and angle < pi * 0.9:
                    return atom.id, atom_id1, atom_id2
    return tuple(tri)


def _normalize_rotation(fragment, tri):
    tri0, tri1, tri2 = tri
    atom_by_id = {atom.id: atom for atom in fragment.atoms}
    moved = fragment.copy()
    mv = atom_by_id[tri0]
    moved.translate(Vector3d(-mv.x, -mv.y, -mv.z))
    moved_by_id = {atom.id: atom for atom in moved.atoms}
    theta = phi = psi = 0.0
    if tri1 == -1:
        return theta, phi, psi
    if tri2 == -1:
        vec = moved_by_id[tri1]
        phi = -_angle(Vector3d(0, 1, 0), Vector3d(0, vec.y, vec.z))
        theta = -_angle(Vector3d(1, 0, 0), Vector3d(vec.x, math.sqrt(vec.y * vec.y + vec.z * vec.z), 0))
        if vec.z < 0:
            phi = -phi
    else:
        vec1 = moved_by_id[tri1]
        vec2 = vec1.cross(moved_by_id[tri2])
        psi = _angle(Vector3d(0, 1, 0), Vector3d(vec2.x, vec2.y, 0))
        phi = _angle(Vector3d(0, 0, 1), Vector3d(0, math.sqrt(vec2.x * vec2.x + vec2.y * vec2.y), vec2.z))
        if vec2.x < 0:
            psi = -psi
        rotated_vec1 = vec1.rotated(0, phi, psi)
        theta = -_angle(Vector3d(1, 0, 0), rotated_vec1)
        if rotated_vec1.y < 0:
            theta = -theta
    return theta, phi, psi


def _normalized_fragment(fragment):
    normalized = fragment.copy()
    center = normalized.center()
    normalized.translate(Vector3d(-center.x, -center.y, -center.z))
    tri = _fragment_reference_triplet(normalized)
    theta, phi, psi = _normalize_rotation(normalized, tri)
    normalized.rotate(theta, phi, psi)
    return normalized


def _fragment_grid_storage_size(conf, score_num):
    count = int((conf.mem_size * 1024 * 1024) / (int(score_num.x) * int(score_num.y) * int(score_num.z) * 4))
    if conf.reuse_grid == ReuseStrategy.NONE:
        return 1
    return max(1, count)


def _prepare_fragment_library(ligands, reorder=True):
    fragvecs = []
    frag_library = []
    frag_importance = []
    fragmap = {}
    for ligand in ligands:
        fragments = decompose_molecule(ligand)
        ligand_fragvecs = []
        for fragment in fragments:
            ob_fragment = _to_obmol(fragment)
            fragment.smiles = _canonical_smiles(ob_fragment)
            labels = canonical_labels(ob_fragment)
            fragment.renumbering(fragment.size, labels)
            signature = fragment.smiles or _fragment_signature(fragment)
            if signature in fragmap:
                frag_idx = fragmap[signature]
                frag_importance[frag_idx] += fragment.size
            else:
                frag_idx = len(frag_library)
                fragmap[signature] = frag_idx
                normalized = _normalized_fragment(fragment)
                normalized.idx = frag_idx
                frag_library.append(normalized)
                frag_importance.append(0)
            fragment.idx = frag_idx
            ligand_fragvecs.append({"pos": fragment.center(), "frag_idx": frag_idx, "size": fragment.size, "rank": 0})
        fragvecs.append(ligand_fragvecs)

    sorted_lig = list(range(len(ligands)))
    if reorder and frag_library:
        fragrank = list(range(len(frag_library)))
        fragrank.sort(key=lambda idx: frag_importance[idx], reverse=True)
        inverse = [0 for _ in fragrank]
        for rank, frag_idx in enumerate(fragrank):
            inverse[frag_idx] = rank
        for ligand_fragvecs in fragvecs:
            for item in ligand_fragvecs:
                item["rank"] = inverse[item["frag_idx"]]
            ligand_fragvecs.sort(key=lambda item: item["rank"])
        sorted_lig.sort(key=lambda lig_idx: [(item["rank"], item["frag_idx"], item["size"]) for item in fragvecs[lig_idx]])
    return frag_library, fragvecs, sorted_lig


def _sample_fragment_grid_for_search(fragment_grid, start, offset, ratio, search_num):
    shape = (int(search_num.x), int(search_num.y), int(search_num.z))
    out = np.full(shape, LIMIT_ENERGY, dtype=np.float32)
    values = fragment_grid.values3d()
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


def _dock_with_fragment_grids(config_path, conf, receptor, ligands):
    atom_grids = _read_atom_grids(config_path, conf)
    score_num = atom_grids[0].num
    search_num = _search_grid_num(conf)
    ratio = _grid_ratio(conf)
    search_grid = InterEnergyGrid(atom_grids[0].center, conf.grid.search_pitch, search_num, 0.0)
    distance_grid = make_distance_grid(atom_grids[0].center, atom_grids[0].pitch, score_num, receptor)
    optimizer = OptimizerGrid(atom_grids, conf.local_max_rmsd)
    ligand_rotations = make_initial_rotations(filename=conf.rotangs_file or None)
    fragment_rotations = make_rotations_60()
    ec = EnergyCalculator(1.0)

    bases = []
    for ligand in ligands:
        base = ligand.copy()
        base.delete_hydrogens()
        base.intra_energy = ec.calc_intra_energy(base)
        center = base.center()
        base.translate(Vector3d(-center.x, -center.y, -center.z))
        bases.append(base)

    frag_library, fragvecs, sorted_lig = _prepare_fragment_library(bases, conf.reorder)
    cache_size = _fragment_grid_storage_size(conf, score_num)
    if conf.reuse_grid == ReuseStrategy.OFFLINE:
        sequences = [[(item["frag_idx"], item["size"]) for item in fragvecs[idx]] for idx in sorted_lig]
        schedule = build_offline_schedule(sequences, cache_size, len(frag_library))
        container = FragmentInterEnergyGridContainer(cache_size, schedule)
    elif conf.reuse_grid == ReuseStrategy.ONLINE:
        container = FragmentInterEnergyGridContainer(cache_size)
    else:
        container = FragmentInterEnergyGridContainer(1)

    candidates_by_ligand = {ligand.identifier: [] for ligand in bases}
    score_start = _to_score_num(0, score_num, search_num, ratio)

    for lig_idx in sorted_lig:
        base = bases[lig_idx]
        frag_list = fragvecs[lig_idx]
        scores = [
            InterEnergyGrid(atom_grids[0].center, conf.grid.search_pitch, search_num, base.intra_energy)
            for _ in ligand_rotations
        ]
        rel_pos = []
        for rotation in ligand_rotations:
            offsets = []
            for item in frag_list:
                rotated_pos = item["pos"].rotated(rotation.x, rotation.y, rotation.z)
                offsets.append(_round_fragment_offset(rotated_pos, conf.grid.score_pitch))
            rel_pos.append(offsets)

        for frag_idx, item in enumerate(frag_list):
            fragid = item["frag_idx"]
            if not container.is_registered(fragid):
                container.insert(FragmentInterEnergyGrid(frag_library[fragid], fragment_rotations, atom_grids, distance_grid))
            fg = container.get(fragid)
            container.next()

            for rotid in range(len(ligand_rotations)):
                offset = rel_pos[rotid][frag_idx]
                scores[rotid].values3d()[:] += _sample_fragment_grid_for_search(fg.grid, score_start, offset, ratio, search_num)

        identifier = base.identifier
        best_param = None
        for rotid, rotation in enumerate(ligand_rotations):
            for x in range(int(search_num.x)):
                for y in range(int(search_num.y)):
                    for z in range(int(search_num.z)):
                        grid_score = scores[rotid].get_inter_energy(x, y, z)
                        if best_param is None or grid_score < best_param[0]:
                            best_param = (grid_score, rotid, x, y, z)
                        if grid_score < conf.output_score_threshold:
                            candidates_by_ligand[identifier].append((grid_score, lig_idx, rotid, x, y, z))
        if not candidates_by_ligand[identifier] and best_param is not None:
            grid_score, rotid, x, y, z = best_param
            candidates_by_ligand[identifier].append((grid_score, lig_idx, rotid, x, y, z))
        if len(candidates_by_ligand[identifier]) > int(conf.poses_per_lig_before_opt) * 2:
            candidates_by_ligand[identifier].sort(key=lambda item: item[0])
            del candidates_by_ligand[identifier][int(conf.poses_per_lig_before_opt) :]

    scored = []
    for candidates in candidates_by_ligand.values():
        candidates.sort(key=lambda item: item[0])
        out_candidates = []
        for _, source_lig_idx, rotid, x, y, z in candidates[: int(conf.poses_per_lig_before_opt)]:
            base = bases[source_lig_idx]
            rotation = ligand_rotations[rotid]
            pose = base.copy()
            pose.rotate(rotation.x, rotation.y, rotation.z)
            pose.translate(search_grid.convert(x, y, z))
            total = optimizer.calc_total_energy(pose)
            if not conf.no_local_opt:
                total = optimizer.optimize(pose)
            inter = total - base.intra_energy
            score = inter / (1.0 + 0.05846 * base.nrots())
            out_candidates.append((score, pose))
        out_candidates.sort(key=lambda item: item[0])
        scored.extend(_select_output_poses(out_candidates, conf.pose_min_rmsd, max(1, int(conf.poses_per_lig))))
    scored.sort(key=lambda item: item[0])
    return [(mol, score) for score, mol in scored]


def _select_output_poses(scored_poses, pose_min_rmsd, max_poses):
    selected = []
    for score, pose in scored_poses:
        if len(selected) >= max_poses:
            break
        min_rmsd = min((pose.calc_rmsd(existing_pose) for _, existing_pose in selected), default=float("inf"))
        if min_rmsd > pose_min_rmsd:
            selected.append((score, pose))
    return selected


def _dock_with_atom_grids(config_path, conf, ligands):
    atom_grids = _read_atom_grids(config_path, conf)
    optimizer = OptimizerGrid(atom_grids, conf.local_max_rmsd)
    search_grid = InterEnergyGrid(atom_grids[0].center, conf.grid.search_pitch, _search_grid_num(conf), 0.0)
    rotations = make_initial_rotations(filename=conf.rotangs_file or None)
    ec = EnergyCalculator(1.0)
    scored = []
    for ligand in ligands:
        base = ligand.copy()
        base.delete_hydrogens()
        base.intra_energy = ec.calc_intra_energy(base)
        if not (conf.score_only or conf.local_only):
            center = base.center()
            base.translate(Vector3d(-center.x, -center.y, -center.z))

        candidates = []
        if conf.score_only or conf.local_only:
            total = _score_on_atom_grids(base, atom_grids) + base.intra_energy
            if conf.local_only and not conf.no_local_opt:
                total = optimizer.optimize(base)
            inter = total - base.intra_energy
            candidates.append((inter / (1.0 + 0.05846 * base.nrots()), base))
        else:
            for rotid, rotation in enumerate(rotations):
                rotated = base.copy()
                rotated.rotate(rotation.x, rotation.y, rotation.z)
                for x in range(int(search_grid.num.x)):
                    for y in range(int(search_grid.num.y)):
                        for z in range(int(search_grid.num.z)):
                            pose = rotated.copy()
                            pose.translate(search_grid.convert(x, y, z))
                            total = _score_on_atom_grids(pose, atom_grids) + base.intra_energy
                            if total < conf.output_score_threshold:
                                if not conf.no_local_opt:
                                    total = optimizer.optimize(pose)
                                inter = total - base.intra_energy
                                score = inter / (1.0 + 0.05846 * base.nrots())
                                candidates.append((score, pose))
            if not candidates:
                best = None
                for rotation in rotations:
                    rotated = base.copy()
                    rotated.rotate(rotation.x, rotation.y, rotation.z)
                    for x in range(int(search_grid.num.x)):
                        for y in range(int(search_grid.num.y)):
                            for z in range(int(search_grid.num.z)):
                                pose = rotated.copy()
                                pose.translate(search_grid.convert(x, y, z))
                                total = _score_on_atom_grids(pose, atom_grids) + base.intra_energy
                                if not conf.no_local_opt:
                                    total = optimizer.optimize(pose)
                                inter = total - base.intra_energy
                                score = inter / (1.0 + 0.05846 * base.nrots())
                                if best is None or score < best[0]:
                                    best = (score, pose)
                if best is not None:
                    candidates.append(best)
        candidates.sort(key=lambda item: item[0])
        scored.extend(candidates[: max(1, int(conf.poses_per_lig))])
    scored.sort(key=lambda item: item[0])
    return [(mol, score) for score, mol in scored]


def build_parser():
    parser = argparse.ArgumentParser(prog="restretto")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func in [
        ("atomgrid-gen", atomgrid_gen),
        ("score-only", score_only),
        ("intraenergy-only", intraenergy_only),
        ("atom-docking", conformer_docking),
        ("easytest-docking", conformer_docking),
    ]:
        cmd = sub.add_parser(name)
        cmd.add_argument("config")
        cmd.add_argument("--grid", "-g")
        cmd.add_argument("--receptor", "-r")
        cmd.add_argument("--ligand", "-l", nargs="+")
        cmd.add_argument("--output", "-o")
        cmd.add_argument("--memsize", "-m", type=int)
        cmd.set_defaults(func=func)
    dock = sub.add_parser("conformer-docking")
    dock.add_argument("--full-rotation", action="store_true")
    dock.add_argument("--output", "-o")
    dock.add_argument("--ligand", "-l", nargs="+")
    dock.add_argument("--receptor", "-r")
    dock.add_argument("--grid", "-g")
    dock.add_argument("--memsize", "-m", type=int)
    dock.add_argument("--score-only", action="store_true")
    dock.add_argument("--no-local-opt", action="store_true")
    dock.add_argument("--local-only", action="store_true")
    dock.add_argument("--local-max-rmsd", type=float)
    dock.add_argument("--log")
    dock.add_argument("--poses-per-lig", type=int)
    dock.add_argument("--min-rmsd", type=float)
    dock.add_argument("--poses-per-lig-before-opt", type=int)
    dock.add_argument("--output-score-threshold", type=float)
    dock.add_argument("--dxgrid")
    dock.add_argument("config")
    dock.set_defaults(func=conformer_docking)
    decomp = sub.add_parser("decompose")
    decomp.add_argument("config", nargs="?")
    decomp.add_argument("--ligand", "-l", nargs="+")
    decomp.add_argument("--fragment", "-f")
    decomp.add_argument("--output", required=True)
    decomp.add_argument("--log")
    decomp.add_argument("--capping_atomic_num", type=int, default=-1)
    decomp.add_argument("--enable_carbon_capping", action="store_true")
    decomp.add_argument("--ins_fragment_id", action="store_true")
    decomp.add_argument("--max_ring_size", type=int, default=-1)
    decomp.add_argument("--no_merge_solitary", action="store_true")
    decomp.set_defaults(func=decompose)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
