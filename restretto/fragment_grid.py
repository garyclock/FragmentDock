from __future__ import annotations

from dataclasses import dataclass
from queue import Queue

import numpy as np

from .constants import LIMIT_ENERGY, XS_TYPE_H
from .geometry import round_half_up
from .grid import InterEnergyGrid


@dataclass
class _MCFPNode:
    to: int
    cost: int


@dataclass
class _FlowNode:
    to: int = -1
    from_: int = -1
    cap: int = 0
    rev_cap: int = 0
    rev_one_cap: int = 0


class _LeftBackSSP:
    INFTY = 1000000007
    CASE_BACKWARD = -1
    CASE_STEPBACKWARD = -2
    CASE_STEPFORWARD = -3
    CASE_UNKNOWN = -4

    def __init__(self, graph):
        self.graph = graph
        self.flow_nodes = [_FlowNode() for _ in graph]
        self.prev = []

    def initialize_flow_nodes(self):
        for i, node in enumerate(self.graph):
            fnode = self.flow_nodes[i]
            fnode.to = node.to
            fnode.rev_cap = 0
            fnode.rev_one_cap = 0
            if fnode.to >= 0:
                fnode.cap = 1
                self.flow_nodes[fnode.to].from_ = i
            else:
                fnode.cap = 0

    def iterative_flow(self):
        if not self.prev:
            self.prev = [self.CASE_UNKNOWN for _ in self.graph]
        dist = [self.INFTY for _ in self.graph]
        dist[0] = 0
        now = 0
        while now < len(self.graph):
            nxt = now + 1
            cost = dist[now]
            if self.flow_nodes[now].cap:
                to = self.flow_nodes[now].to
                to_cost = cost - self.graph[now].cost
                if dist[to] > to_cost:
                    dist[to] = to_cost
                    self.prev[to] = now
            if now + 1 < len(self.graph):
                to = now + 1
                to_cost = cost
                if dist[to] > to_cost:
                    dist[to] = to_cost
                    self.prev[to] = self.CASE_STEPFORWARD
            if self.prev[now] != self.CASE_STEPFORWARD and self.flow_nodes[now].rev_one_cap:
                to = now - 1
                to_cost = cost
                if dist[to] > to_cost:
                    dist[to] = to_cost
                    self.prev[to] = self.CASE_STEPBACKWARD
                    nxt = to
            if self.flow_nodes[now].rev_cap:
                to = self.flow_nodes[now].from_
                to_cost = cost + self.graph[to].cost
                if dist[to] > to_cost:
                    dist[to] = to_cost
                    self.prev[to] = self.CASE_BACKWARD
                    nxt = to
            now = nxt

        now = len(self.graph) - 1
        while now != 0:
            prev = self.prev[now]
            if prev == self.CASE_BACKWARD:
                self.flow_nodes[now].cap += 1
                now = self.flow_nodes[now].to
                self.flow_nodes[now].rev_cap -= 1
            elif prev == self.CASE_STEPBACKWARD:
                now = now + 1
                self.flow_nodes[now].rev_one_cap -= 1
            elif prev == self.CASE_STEPFORWARD:
                self.flow_nodes[now].rev_one_cap += 1
                now = now - 1
            else:
                self.flow_nodes[now].rev_cap += 1
                now = prev
                self.flow_nodes[now].cap -= 1

        return dist[-1]

    def run(self, supply):
        total = 0
        self.initialize_flow_nodes()
        for _ in range(supply):
            total += self.iterative_flow()
        return total


def _run_left_back_ssp(graph, supply):
    if not graph:
        return []
    solver = _LeftBackSSP(graph)
    solver.run(supply - 1)
    ret = [-1 for _ in graph]
    empty = Queue()
    for i in range(supply):
        empty.put(i)
    for i in range(len(graph)):
        if solver.flow_nodes[i].to == i:
            ret[i] = ret[i - 1]
        elif solver.flow_nodes[i].rev_cap:
            ret[i] = ret[solver.flow_nodes[i].from_ - 1]
        else:
            ret[i] = empty.get()
        if i == len(graph) - 1 or (
            solver.flow_nodes[i + 1].to != i + 1
            and (solver.flow_nodes[i + 1].cap == 1 or solver.flow_nodes[i + 1].to < 0)
        ):
            empty.put(ret[i])
    return ret


def build_offline_schedule(fragment_sequences, cache_size, fragment_count):
    """Port of fraggrid_main.cc makeGraph plus CalcMCFP::runLeftBackSSP."""
    cache_size = int(cache_size)
    if cache_size <= 0:
        raise ValueError("cache_size must be positive")
    size = sum(len(sequence) for sequence in fragment_sequences)
    before = [-1 for _ in range(fragment_count)]
    graph = [_MCFPNode(-1, 0) for _ in range(size)]
    cursor = 0
    for sequence in fragment_sequences:
        for frag_idx, frag_size in sequence:
            if before[frag_idx] == cursor:
                graph[cursor] = _MCFPNode(cursor, 0)
            elif before[frag_idx] != -1:
                graph[before[frag_idx]] = _MCFPNode(cursor, frag_size)
            cursor += 1
            before[frag_idx] = cursor
    return _run_left_back_ssp(graph, cache_size)


def make_distance_grid(center, pitch, num, receptor_mol):
    grid = InterEnergyGrid(center, pitch, num)
    shape = (int(grid.num.x), int(grid.num.y), int(grid.num.z))
    xs = (np.arange(shape[0], dtype=np.float32) - (shape[0] - 1) / 2.0) * pitch.x + center.x
    ys = (np.arange(shape[1], dtype=np.float32) - (shape[1] - 1) / 2.0) * pitch.y + center.y
    zs = (np.arange(shape[2], dtype=np.float32) - (shape[2] - 1) / 2.0) * pitch.z + center.z
    xvals, yvals, zvals = np.meshgrid(xs, ys, zs, indexing="ij")
    min_dist = np.full(shape, 1e4, dtype=np.float32)
    for atom in receptor_mol.atoms:
        if atom.xs_type == XS_TYPE_H:
            continue
        dist = np.sqrt((xvals - atom.x) ** 2 + (yvals - atom.y) ** 2 + (zvals - atom.z) ** 2, dtype=np.float32)
        np.minimum(min_dist, dist, out=min_dist)
    grid.grid = min_dist.reshape(-1).astype(np.float32, copy=False)
    return grid


def _shifted_grid_values(grid, dx, dy, dz, shape):
    values = grid.values3d()
    out = np.full(shape, LIMIT_ENERGY, dtype=np.float32)
    src_x0 = max(0, dx)
    src_y0 = max(0, dy)
    src_z0 = max(0, dz)
    src_x1 = min(shape[0], shape[0] + dx)
    src_y1 = min(shape[1], shape[1] + dy)
    src_z1 = min(shape[2], shape[2] + dz)
    if src_x0 >= src_x1 or src_y0 >= src_y1 or src_z0 >= src_z1:
        return out
    dst_x0 = src_x0 - dx
    dst_y0 = src_y0 - dy
    dst_z0 = src_z0 - dz
    dst_x1 = src_x1 - dx
    dst_y1 = src_y1 - dy
    dst_z1 = src_z1 - dz
    out[dst_x0:dst_x1, dst_y0:dst_y1, dst_z0:dst_z1] = values[src_x0:src_x1, src_y0:src_y1, src_z0:src_z1]
    return out


class FragmentInterEnergyGrid:
    def __init__(self, orig_frag, rot_angles, atom_grids, distance_grid):
        if not atom_grids:
            raise ValueError("atom_grids is empty")
        self.frag_idx = getattr(orig_frag, "idx", -1)
        num = atom_grids[0].num
        self.grid = InterEnergyGrid(atom_grids[0].center, atom_grids[0].pitch, num, LIMIT_ENERGY)
        rotations = list(rot_angles)
        if orig_frag.size <= 1:
            rotations = rotations[:1] or []
        radius = orig_frag.radius()
        shape = (int(num.x), int(num.y), int(num.z))
        best_values = self.grid.values3d()
        distance_values = distance_grid.values3d()

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
            np.minimum(best_values, total, out=best_values)


class FragmentInterEnergyGridContainer:
    def __init__(self, size=0, indices_to_save=None):
        self.size = int(size)
        self.grids = [None for _ in range(self.size)]
        self.indices_to_save = list(indices_to_save or [])
        self.step = 0
        self.last_used = [-1 for _ in range(self.size)]
        self.offline = indices_to_save is not None
        if self.offline:
            for index in self.indices_to_save:
                if index >= self.size:
                    raise ValueError("Invalid index to save")

    def _search(self, fragid):
        for idx, grid in enumerate(self.grids):
            if grid is not None and grid.frag_idx == fragid:
                return idx
        return -1

    def _lru_idx(self):
        if not self.last_used:
            raise ValueError("empty fragment grid container")
        return min(range(self.size), key=lambda idx: self.last_used[idx])

    def insert(self, grid):
        if self.offline:
            self.grids[self.indices_to_save[self.step]] = grid
        else:
            idx = self._lru_idx()
            self.grids[idx] = grid
            self.last_used[idx] = self.step

    def is_registered(self, fragid):
        if self.offline:
            grid = self.grids[self.indices_to_save[self.step]]
            return grid is not None and grid.frag_idx == fragid
        return self._search(fragid) != -1

    def get(self, fragid):
        if not self.is_registered(fragid):
            raise KeyError(fragid)
        if self.offline:
            return self.grids[self.indices_to_save[self.step]]
        idx = self._search(fragid)
        self.last_used[idx] = self.step
        return self.grids[idx]

    def next(self):
        self.step += 1
