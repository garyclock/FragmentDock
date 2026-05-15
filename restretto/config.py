from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from .constants import EPS
from .geometry import Point3d, ceili


class ReuseStrategy(Enum):
    OFFLINE = "offline"
    ONLINE = "online"
    NONE = "none"


@dataclass
class SearchGrid:
    center: Point3d = field(default_factory=lambda: Point3d(0.0, 0.0, 0.0))
    outer_width: Point3d = field(default_factory=lambda: Point3d(0.0, 0.0, 0.0))
    inner_width: Point3d = field(default_factory=lambda: Point3d(0.0, 0.0, 0.0))
    search_pitch: Point3d = field(default_factory=lambda: Point3d(1.0, 1.0, 1.0))
    score_pitch: Point3d = field(default_factory=lambda: Point3d(1.0, 1.0, 1.0))


@dataclass
class DockingConfiguration:
    grid: SearchGrid = field(default_factory=SearchGrid)
    ligand_files: list = field(default_factory=list)
    receptor_file: str = ""
    output_file: str = "out.sdf"
    log_file: str = ""
    grid_folder: str = ""
    dxgrid_folder: str = ""
    rotangs_file: str = ""
    reuse_grid: ReuseStrategy = ReuseStrategy.OFFLINE
    reorder: bool = True
    mem_size: int = 8000
    poses_per_lig: int = 1
    poses_per_lig_before_opt: int = 2000
    output_score_threshold: float = -3.0
    pose_min_rmsd: float = 0.5
    no_local_opt: bool = False
    score_only: bool = False
    local_only: bool = False
    local_max_rmsd: float = 1e10
    rad_scale: float = 0.95

    def check_validity(self):
        for score, search in zip(self.grid.score_pitch.as_tuple(), self.grid.search_pitch.as_tuple()):
            ratio = round(search / score)
            if abs(score * ratio - search) >= EPS:
                raise ValueError("SEARCH_PITCH / SCORING_PITCH must be an integer ratio")
        num = tuple(ceili(width / 2.0 / pitch) * 2 + 1 for width, pitch in zip(self.grid.outer_width.as_tuple(), self.grid.score_pitch.as_tuple()))
        grid_count = int((self.mem_size * 1024 * 1024) / (num[0] * num[1] * num[2] * 4))
        if grid_count <= 0:
            raise ValueError("MEMORY_SIZE is too small to store one fragment grid")
        if self.score_only and self.local_only:
            raise ValueError("SCORE_ONLY and LOCAL_ONLY cannot both be true")


def _parse_point(value):
    parts = [float(part.strip()) for part in value.split(",")]
    if len(parts) != 3:
        raise ValueError("expected three comma-separated values")
    return Point3d(parts[0], parts[1], parts[2])


def _parse_bool(value):
    lowered = value.strip().lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    raise ValueError("expected true or false")


def parse_in_file(filename):
    conf = DockingConfiguration()
    for raw in Path(filename).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        key, _, value = line.partition(" ")
        value = value.strip()
        if key == "OUTERBOX":
            conf.grid.outer_width = _parse_point(value)
        elif key == "INNERBOX":
            conf.grid.inner_width = _parse_point(value)
        elif key == "BOX_CENTER":
            conf.grid.center = _parse_point(value)
        elif key == "SEARCH_PITCH":
            conf.grid.search_pitch = _parse_point(value)
        elif key == "SCORING_PITCH":
            conf.grid.score_pitch = _parse_point(value)
        elif key == "REUSE_FRAG_GRID":
            lowered = value.lower()
            if lowered in {"false", "none", "noreuse"}:
                conf.reuse_grid = ReuseStrategy.NONE
            elif lowered == "online":
                conf.reuse_grid = ReuseStrategy.ONLINE
            else:
                conf.reuse_grid = ReuseStrategy.OFFLINE
        elif key == "REORDER_LIGANDS":
            conf.reorder = _parse_bool(value)
        elif key == "MEMORY_SIZE":
            conf.mem_size = int(value)
        elif key == "RECEPTOR":
            conf.receptor_file = value
        elif key == "LIGAND":
            conf.ligand_files.append(value)
        elif key == "OUTPUT":
            conf.output_file = value
        elif key == "LOG":
            conf.log_file = value
        elif key == "GRID_FOLDER":
            conf.grid_folder = value
        elif key == "DXGRID_FOLDER":
            conf.dxgrid_folder = value
        elif key == "ROTANGS":
            conf.rotangs_file = value
        elif key == "POSES_PER_LIG":
            conf.poses_per_lig = int(value)
        elif key == "POSES_PER_LIG_BEFORE_OPT":
            conf.poses_per_lig_before_opt = int(value)
        elif key == "OUTPUT_SCORE_THRESHOLD":
            conf.output_score_threshold = float(value)
        elif key == "MIN_RMSD":
            conf.pose_min_rmsd = float(value)
        elif key == "NO_LOCAL_OPT":
            conf.no_local_opt = _parse_bool(value)
        elif key == "SCORE_ONLY":
            conf.score_only = _parse_bool(value)
        elif key == "LOCAL_ONLY":
            conf.local_only = _parse_bool(value)
        elif key == "LOCAL_MAX_RMSD":
            conf.local_max_rmsd = float(value)
    conf.check_validity()
    return conf
