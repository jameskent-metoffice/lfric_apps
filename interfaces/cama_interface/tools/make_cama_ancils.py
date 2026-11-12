#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Create the ancillary files needed to couple CaMa-Flood to LFRic and NEMO.

Sub-commands:

river-numbers
    Create the map of NEMO river numbers on the CaMa-Flood grid. Every
    CaMa-Flood river mouth that drains to the ocean is given the NEMO river
    number found at its location in the JULES river number ancillary, which is
    also what NEMO uses to receive its runoff. Where the JULES ancillary has no
    river at that location, the nearest river within a search radius is used.
    The file holds one 4 byte integer per CaMa-Flood grid point, x varying
    fastest with the first row in the north, like the CaMa-Flood binary maps.
    Inland basins are given 0, so their water stays in CaMa-Flood.

diminfo
    Write the CaMa-Flood dimension file (CDIMINFO) for runoff supplied on the
    CaMa-Flood map grid itself (no input matrix).

scrip-grid
    Write a SCRIP description of the CaMa-Flood grid, with rows ordered from
    north to south as in OASIS, for generating the LFRic to CaMa-Flood
    remapping weights.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

CAMA_OCEAN_MOUTH = -9


def read_nextx(map_dir, nx, ny, big_endian=False):
    """Read the downstream x index of every CaMa-Flood grid point."""
    dtype = ">i4" if big_endian else "<i4"
    path = Path(map_dir) / "nextxy.bin"
    data = np.fromfile(path, dtype=dtype)
    if data.size != 2 * nx * ny:
        raise ValueError(
            f"{path} holds {data.size} values, expected {2 * nx * ny} "
            f"for a {nx}x{ny} grid"
        )
    return data[: nx * ny].reshape((ny, nx))


def cell_centres(west, east, north, south, nx, ny):
    """Longitudes and latitudes of the cell centres, first row in the north."""
    dlon = (east - west) / nx
    dlat = (north - south) / ny
    lon = west + (np.arange(nx) + 0.5) * dlon
    lat = north - (np.arange(ny) + 0.5) * dlat
    return lon, lat


def lookup_river_numbers(
    mouth_lon, mouth_lat, rivers, riv_lon, riv_lat, search_radius
):
    """
    Find the river number at each mouth location.

    rivers is a 2D array (lat, lon) of river numbers, 0 meaning no river, on
    a regular grid with 1D, monotonic cell centre coordinates riv_lon and
    riv_lat in degrees. A mouth in a cell without a river takes the nearest
    river within search_radius degrees.
    """
    result = np.zeros(len(mouth_lon), dtype=np.int32)
    lon_step = np.abs(np.mean(np.diff(riv_lon)))
    lat_step = np.abs(np.mean(np.diff(riv_lat)))
    river_j, river_i = np.nonzero(rivers > 0)
    river_lon = riv_lon[river_i]
    river_lat = riv_lat[river_j]

    for n, (lon, lat) in enumerate(zip(mouth_lon, mouth_lat)):
        i = int(np.argmin(np.abs(riv_lon - lon)))
        j = int(np.argmin(np.abs(riv_lat - lat)))
        if rivers[j, i] > 0 and (
            abs(riv_lon[i] - lon) <= lon_step and abs(riv_lat[j] - lat) <= lat_step
        ):
            result[n] = rivers[j, i]
            continue
        if river_lon.size == 0:
            continue
        # Distance in degrees, with longitude scaled by latitude
        dlon = (river_lon - lon + 180.0) % 360.0 - 180.0
        dlat = river_lat - lat
        dist = np.hypot(dlon * np.cos(np.radians(lat)), dlat)
        nearest = int(np.argmin(dist))
        if dist[nearest] <= search_radius:
            result[n] = rivers[river_j[nearest], river_i[nearest]]
    return result


def read_jules_river_numbers(path, variable, lon_name, lat_name):
    """Read the JULES river number ancillary as (lat, lon) with 1D coordinates."""
    from netCDF4 import Dataset  # pylint: disable=import-outside-toplevel

    with Dataset(path) as nc:
        var = nc.variables[variable]
        data = np.ma.filled(np.squeeze(var[:]), 0).astype(np.int64)
        if data.ndim != 2:
            raise ValueError(f"{variable} in {path} is not a 2D field")
        lon = np.asarray(nc.variables[lon_name][:], dtype=float)
        lat = np.asarray(nc.variables[lat_name][:], dtype=float)
    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError(
            "The river number file needs 1D longitude and latitude "
            f"coordinates, found {lon_name} and {lat_name} with "
            f"{lon.ndim} and {lat.ndim} dimensions"
        )
    if data.shape != (lat.size, lon.size):
        raise ValueError(
            f"{variable} has shape {data.shape}, expected {(lat.size, lon.size)}"
        )
    return data, lon, lat


def river_numbers_command(args):
    """Create the river number map on the CaMa-Flood grid."""
    nextx = read_nextx(args.map_dir, args.nx, args.ny, args.big_endian)
    lon, lat = cell_centres(
        args.west, args.east, args.north, args.south, args.nx, args.ny
    )
    mouth_j, mouth_i = np.nonzero(nextx == CAMA_OCEAN_MOUTH)

    rivers, riv_lon, riv_lat = read_jules_river_numbers(
        args.river_number_file, args.variable, args.lon_name, args.lat_name
    )
    found = lookup_river_numbers(
        lon[mouth_i], lat[mouth_j], rivers, riv_lon, riv_lat, args.search_radius
    )

    result = np.zeros((args.ny, args.nx), dtype="<i4")
    result[mouth_j, mouth_i] = found
    result.tofile(args.output)

    n_rivers = int(rivers.max())
    represented = np.unique(found[found > 0])
    print(f"{mouth_j.size} ocean river mouths in the CaMa-Flood map")
    print(f"{int((found > 0).sum())} mapped to a river number")
    print(f"{int((found == 0).sum())} not mapped (no river within "
          f"{args.search_radius} degrees); their water is not sent to the ocean")
    print(f"{represented.size} of {n_rivers} river numbers receive a mouth")
    print(f"Written to {args.output}")


def diminfo_command(args):
    """Write the CDIMINFO file for runoff on the CaMa-Flood grid."""
    lines = [
        f"{args.nx:>12}     !! nXX",
        f"{args.ny:>12}     !! nYY",
        f"{args.nlfp:>12}     !! floodplain layer",
        f"{args.nx:>12}     !! input nXX",
        f"{args.ny:>12}     !! input nYY",
        f"{1:>12}     !! input num",
        "NONE",
        f"{args.west:>12.3f}     !! west  edge",
        f"{args.east:>12.3f}     !! east  edge",
        f"{args.north:>12.3f}     !! north edge",
        f"{args.south:>12.3f}     !! south edge",
    ]
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"Written to {args.output}")


def scrip_grid_command(args):
    """Write a SCRIP grid description of the CaMa-Flood grid."""
    from netCDF4 import Dataset  # pylint: disable=import-outside-toplevel

    lon, lat = cell_centres(
        args.west, args.east, args.north, args.south, args.nx, args.ny
    )
    dlon = (args.east - args.west) / args.nx
    dlat = (args.north - args.south) / args.ny
    lon2d, lat2d = np.meshgrid(lon, lat)

    # Corners in counter-clockwise order: south-west, south-east, north-east,
    # north-west
    corner_lon = np.stack(
        [lon2d - dlon / 2, lon2d + dlon / 2, lon2d + dlon / 2, lon2d - dlon / 2],
        axis=-1,
    )
    corner_lat = np.stack(
        [lat2d - dlat / 2, lat2d - dlat / 2, lat2d + dlat / 2, lat2d + dlat / 2],
        axis=-1,
    )

    with Dataset(args.output, "w") as nc:
        nc.createDimension("grid_size", args.nx * args.ny)
        nc.createDimension("grid_corners", 4)
        nc.createDimension("grid_rank", 2)
        dims = nc.createVariable("grid_dims", "i4", ("grid_rank",))
        dims[:] = [args.nx, args.ny]
        for name, values in (("grid_center_lat", lat2d), ("grid_center_lon", lon2d)):
            var = nc.createVariable(name, "f8", ("grid_size",))
            var.units = "degrees"
            var[:] = values.reshape(-1)
        for name, values in (
            ("grid_corner_lat", corner_lat),
            ("grid_corner_lon", corner_lon),
        ):
            var = nc.createVariable(name, "f8", ("grid_size", "grid_corners"))
            var.units = "degrees"
            var[:] = values.reshape(-1, 4)
        mask = nc.createVariable("grid_imask", "i4", ("grid_size",))
        mask[:] = 1
    print(f"Written to {args.output}")


def add_grid_arguments(parser):
    """Arguments describing the CaMa-Flood grid."""
    parser.add_argument("--nx", type=int, default=1440, help="Grid points in x")
    parser.add_argument("--ny", type=int, default=720, help="Grid points in y")
    parser.add_argument("--west", type=float, default=-180.0)
    parser.add_argument("--east", type=float, default=180.0)
    parser.add_argument("--north", type=float, default=90.0)
    parser.add_argument("--south", type=float, default=-90.0)


def parse_args(argv):
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    rivers = sub.add_parser("river-numbers", help="NEMO river numbers on the grid")
    add_grid_arguments(rivers)
    rivers.add_argument("--map-dir", required=True, help="Directory with nextxy.bin")
    rivers.add_argument("--big-endian", action="store_true")
    rivers.add_argument(
        "--river-number-file", required=True, help="JULES river number NetCDF file"
    )
    rivers.add_argument("--variable", default="river_number")
    rivers.add_argument("--lon-name", default="longitude")
    rivers.add_argument("--lat-name", default="latitude")
    rivers.add_argument(
        "--search-radius",
        type=float,
        default=3.0,
        help="Search radius in degrees for a mouth in a cell without a river",
    )
    rivers.add_argument("--output", required=True)
    rivers.set_defaults(func=river_numbers_command)

    diminfo = sub.add_parser("diminfo", help="CaMa-Flood CDIMINFO file")
    add_grid_arguments(diminfo)
    diminfo.add_argument("--nlfp", type=int, default=10, help="Floodplain layers")
    diminfo.add_argument("--output", required=True)
    diminfo.set_defaults(func=diminfo_command)

    scrip = sub.add_parser("scrip-grid", help="SCRIP grid of the CaMa-Flood grid")
    add_grid_arguments(scrip)
    scrip.add_argument("--output", required=True)
    scrip.set_defaults(func=scrip_grid_command)

    return parser.parse_args(argv)


def main(argv=None):
    """Entry point."""
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args.func(args)


if __name__ == "__main__":
    main()
