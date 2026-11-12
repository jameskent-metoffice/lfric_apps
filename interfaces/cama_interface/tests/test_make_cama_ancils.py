##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Tests for make_cama_ancils.py, using small synthetic grids."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import make_cama_ancils as tool  # noqa: E402


def write_nextxy(path, nextx):
    """Write nextxy.bin holding nextx then nexty (nexty is not used)."""
    nexty = np.ones_like(nextx)
    np.concatenate([nextx.reshape(-1), nexty.reshape(-1)]).astype("<i4").tofile(
        Path(path) / "nextxy.bin"
    )


def write_rivers(path, rivers, lon, lat):
    """Write a river number file on a regular grid."""
    with Dataset(path, "w") as nc:
        nc.createDimension("latitude", lat.size)
        nc.createDimension("longitude", lon.size)
        nc.createVariable("latitude", "f8", ("latitude",))[:] = lat
        nc.createVariable("longitude", "f8", ("longitude",))[:] = lon
        nc.createVariable("river_number", "i4", ("latitude", "longitude"))[:] = rivers


class TestLookup(unittest.TestCase):
    """The lookup of river numbers at the river mouths."""

    def setUp(self):
        self.lon = np.array([0.5, 1.5, 2.5, 3.5])
        self.lat = np.array([3.5, 2.5, 1.5, 0.5])
        self.rivers = np.zeros((4, 4), dtype=np.int64)
        self.rivers[0, 0] = 7
        self.rivers[3, 3] = 9

    def test_mouth_inside_a_river_cell(self):
        found = tool.lookup_river_numbers(
            np.array([0.4]), np.array([3.6]), self.rivers, self.lon, self.lat, 3.0
        )
        self.assertEqual(found[0], 7)

    def test_nearest_river_is_used_without_a_river_in_the_cell(self):
        found = tool.lookup_river_numbers(
            np.array([2.5]), np.array([0.5]), self.rivers, self.lon, self.lat, 3.0
        )
        self.assertEqual(found[0], 9)

    def test_no_river_within_the_search_radius(self):
        found = tool.lookup_river_numbers(
            np.array([2.0]), np.array([2.0]), self.rivers, self.lon, self.lat, 0.5
        )
        self.assertEqual(found[0], 0)


class TestCommands(unittest.TestCase):
    """The command line sub-commands."""

    def test_river_numbers_only_maps_ocean_mouths(self):
        nx, ny = 4, 2
        # Row 1 is the north. Cell (1, 1) is an ocean mouth, (3, 1) an inland
        # basin and (2, 2) an ocean mouth far from any river.
        nextx = np.full((ny, nx), 2, dtype=np.int64)
        nextx[0, 0] = tool.CAMA_OCEAN_MOUTH
        nextx[0, 2] = -10
        nextx[1, 1] = tool.CAMA_OCEAN_MOUTH
        with tempfile.TemporaryDirectory() as tmp:
            write_nextxy(tmp, nextx)
            rivers = np.zeros((2, 4), dtype=np.int64)
            rivers[0, 0] = 3
            write_rivers(
                Path(tmp) / "rivers.nc",
                rivers,
                np.array([0.5, 1.5, 2.5, 3.5]),
                np.array([1.5, 0.5]),
            )
            out = Path(tmp) / "numbers.bin"
            tool.main(
                [
                    "river-numbers", "--nx", str(nx), "--ny", str(ny),
                    "--west", "0", "--east", "4", "--north", "2", "--south", "0",
                    "--map-dir", tmp,
                    "--river-number-file", str(Path(tmp) / "rivers.nc"),
                    "--search-radius", "0.5", "--output", str(out),
                ]
            )
            result = np.fromfile(out, dtype="<i4").reshape(ny, nx)
        expected = np.zeros((ny, nx), dtype=np.int32)
        expected[0, 0] = 3
        np.testing.assert_array_equal(result, expected)

    def test_wrong_nextxy_size_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_nextxy(tmp, np.ones((2, 2), dtype=np.int64))
            with self.assertRaises(ValueError):
                tool.read_nextx(tmp, 4, 4)

    def test_diminfo(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "diminfo.txt"
            tool.main(["diminfo", "--output", str(out)])
            values = [line.split()[0] for line in out.read_text().splitlines()]
        self.assertEqual(values[:6], ["1440", "720", "10", "1440", "720", "1"])

    def test_scrip_grid_is_ordered_north_to_south(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "grid.nc"
            tool.main(
                ["scrip-grid", "--nx", "4", "--ny", "2", "--west", "0",
                 "--east", "4", "--north", "2", "--south", "0",
                 "--output", str(out)]
            )
            with Dataset(out) as nc:
                lat = nc.variables["grid_center_lat"][:]
                lon = nc.variables["grid_center_lon"][:]
                dims = nc.variables["grid_dims"][:]
        np.testing.assert_array_equal(dims, [4, 2])
        np.testing.assert_allclose(lon[:4], [0.5, 1.5, 2.5, 3.5])
        np.testing.assert_allclose(lat, [1.5] * 4 + [0.5] * 4)


if __name__ == "__main__":
    unittest.main()
