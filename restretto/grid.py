from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from .constants import INF_ENERGY, LIMIT_ENERGY
from .geometry import Point3d, Vector3d, round_half_up


class InterEnergyGrid:
    def __init__(self, center=None, pitch=None, num=None, initial=INF_ENERGY):
        self.center = center or Point3d(0.0, 0.0, 0.0)
        self.pitch = pitch or Point3d(1.0, 1.0, 1.0)
        self.num = num or Point3d(0, 0, 0)
        self.grid = np.full(int(self.num.x) * int(self.num.y) * int(self.num.z), float(initial), dtype=np.float32)

    def _offset(self, x, y, z):
        return (int(x) * int(self.num.y) + int(y)) * int(self.num.z) + int(z)

    def set_inter_energy(self, x, y, z, value):
        self.grid[self._offset(x, y, z)] = float(value)

    def add_energy(self, x, y, z, value):
        self.grid[self._offset(x, y, z)] += float(value)

    def get_inter_energy(self, x, y=None, z=None):
        if isinstance(x, Vector3d) and y is None and z is None:
            return self.get_inter_energy(self.convert_x(x), self.convert_y(x), self.convert_z(x))
        if x < 0 or y < 0 or z < 0 or x >= self.num.x or y >= self.num.y or z >= self.num.z:
            return LIMIT_ENERGY
        return self.grid[self._offset(x, y, z)]

    def convert_x(self, vec):
        return round_half_up((vec.x - self.center.x) / self.pitch.x + (self.num.x - 1) / 2.0)

    def convert_y(self, vec):
        return round_half_up((vec.y - self.center.y) / self.pitch.y + (self.num.y - 1) / 2.0)

    def convert_z(self, vec):
        return round_half_up((vec.z - self.center.z) / self.pitch.z + (self.num.z - 1) / 2.0)

    def convert(self, x, y, z):
        return Vector3d(
            (x - (self.num.x - 1) / 2.0) * self.pitch.x + self.center.x,
            (y - (self.num.y - 1) / 2.0) * self.pitch.y + self.center.y,
            (z - (self.num.z - 1) / 2.0) * self.pitch.z + self.center.z,
        )

    def values3d(self):
        return self.grid.reshape((int(self.num.x), int(self.num.y), int(self.num.z)))

    def write_grid(self, filename):
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            stream.write(struct.pack("<3f", self.center.x, self.center.y, self.center.z))
            stream.write(struct.pack("<3f", self.pitch.x, self.pitch.y, self.pitch.z))
            stream.write(struct.pack("<3i", int(self.num.x), int(self.num.y), int(self.num.z)))
            for value in self.grid:
                stream.write(struct.pack("<f", float(value)))

    @classmethod
    def parse_grid(cls, filename):
        with Path(filename).open("rb") as stream:
            center = Point3d(*struct.unpack("<3f", stream.read(12)))
            pitch = Point3d(*struct.unpack("<3f", stream.read(12)))
            num = Point3d(*struct.unpack("<3i", stream.read(12)))
            grid = cls(center, pitch, num)
            count = int(num.x) * int(num.y) * int(num.z)
            grid.grid = np.frombuffer(stream.read(4 * count), dtype="<f4").astype(np.float32, copy=True)
            return grid
