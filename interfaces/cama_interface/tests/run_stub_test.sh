#!/usr/bin/env bash
##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# End to end test of the CaMa-Flood OASIS driver without OASIS3-MCT.
#
# Builds cama.exe against stand-in OASIS and MPI modules (tests/stubs), and runs
# it on a small CaMa-Flood river map. The stub OASIS supplies a uniform runoff
# and records what the driver sends, so the test checks the driver, the CaMa
# namelist set-up and the mapping of river mouths onto river numbers.
#
# Usage:
#   CAMA_SOURCE=<CaMa-Flood_v4 checkout> CAMA_MAP_DIR=<map directory> \
#       run_stub_test.sh
#
# CAMA_MAP_DIR must contain nextxy.bin, ctmare.bin, elevtn.bin, nxtdst.bin,
# rivlen.bin, fldhgt.bin, rivwth.bin, rivhgt.bin, rivman.bin and bifprm.txt,
# for example map/mekong_15min in the CaMa-Flood repository. The grid size
# and bounds are given by CAMA_NX, CAMA_NY, CAMA_WEST, CAMA_EAST, CAMA_NORTH
# and CAMA_SOUTH (the defaults describe the Mekong example).
#
# Requires gfortran, python3 with numpy and netCDF4.
##############################################################################
set -euo pipefail

: "${CAMA_SOURCE:?CAMA_SOURCE must point at the CaMa-Flood_v4 source}"
: "${CAMA_MAP_DIR:?CAMA_MAP_DIR must point at a CaMa-Flood map directory}"
nx=${CAMA_NX:-80}
ny=${CAMA_NY:-120}
west=${CAMA_WEST:-90}
east=${CAMA_EAST:-110}
north=${CAMA_NORTH:-35}
south=${CAMA_SOUTH:-5}
n_rivers=5
runoff=1.0e-5      # kg m-2 s-1
hours=${CAMA_TEST_HOURS:-12}

test_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
interface_dir=$(cd "$test_dir/.." && pwd)
work_dir=${CAMA_TEST_DIR:-$(mktemp -d)}
mkdir -p "$work_dir/stubs" "$work_dir/run"
echo "[INFO] Working in $work_dir"

# Build the driver against the stand-in modules
(cd "$work_dir/stubs" && gfortran -c "$test_dir/stubs/mpi.F90" \
    "$test_dir/stubs/mod_oasis.F90")
CAMA_MAKE_DIR="$work_dir/make" FC=gfortran CAMA_COMPILER=gnu \
    OASIS_INC="$work_dir/stubs" OASIS_LIBS="" \
    OASIS_EXTRA_LIBS="$work_dir/stubs/mpi.o $work_dir/stubs/mod_oasis.o" \
    "$interface_dir/build/build_cama.sh" > "$work_dir/build.log"

cd "$work_dir/run"

# A river number file in which every river drains to river 1..n_rivers
python3 - << EOF
import numpy as np
from netCDF4 import Dataset
lon = np.arange(-179.5, 180, 1.0)
lat = np.arange(-89.5, 90, 1.0)
rivers = (np.arange(lat.size * lon.size).reshape(lat.size, lon.size)
          % $n_rivers + 1).astype("i4")
with Dataset("rivers.nc", "w") as nc:
    nc.createDimension("latitude", lat.size)
    nc.createDimension("longitude", lon.size)
    nc.createVariable("latitude", "f8", ("latitude",))[:] = lat
    nc.createVariable("longitude", "f8", ("longitude",))[:] = lon
    nc.createVariable("river_number", "i4", ("latitude", "longitude"))[:] = rivers
EOF
python3 "$interface_dir/tools/make_cama_ancils.py" river-numbers \
    --nx "$nx" --ny "$ny" --west "$west" --east "$east" --north "$north" \
    --south "$south" --map-dir "$CAMA_MAP_DIR" \
    --river-number-file rivers.nc --output cama_river_number.bin
python3 "$interface_dir/tools/make_cama_ancils.py" diminfo \
    --nx "$nx" --ny "$ny" --west "$west" --east "$east" --north "$north" \
    --south "$south" --output cama_diminfo.txt

cat > cama_oasis.nml << EOF
&cama_oasis_nml
river_number_file='cama_river_number.bin'
cpl_freq=3600
river_count=$n_rivers
water_density=1000.0
/
EOF

cat > input_cmf.nam << EOF
&NRUNVER
LADPSTP=.TRUE.
LPTHOUT=.TRUE.
LROSPLIT=.TRUE.
LRESTART=.FALSE.
/
&NDIMTIME
CDIMINFO='cama_diminfo.txt'
DT=3600
IFRQ_INP=1
/
&NPARAM
/
&NSIMTIME
SYEAR=2000
SMON=1
SDAY=1
SHOUR=0
EYEAR=2000
EMON=1
EDAY=$((1 + hours / 24))
EHOUR=$((hours % 24))
/
&NMAP
LMAPCDF=.FALSE.
CNEXTXY='$CAMA_MAP_DIR/nextxy.bin'
CGRAREA='$CAMA_MAP_DIR/ctmare.bin'
CELEVTN='$CAMA_MAP_DIR/elevtn.bin'
CNXTDST='$CAMA_MAP_DIR/nxtdst.bin'
CRIVLEN='$CAMA_MAP_DIR/rivlen.bin'
CFLDHGT='$CAMA_MAP_DIR/fldhgt.bin'
CRIVWTH='$CAMA_MAP_DIR/rivwth.bin'
CRIVHGT='$CAMA_MAP_DIR/rivhgt.bin'
CRIVMAN='$CAMA_MAP_DIR/rivman.bin'
CPTHOUT='$CAMA_MAP_DIR/bifprm.txt'
/
&NRESTART
CRESTSTO=''
CRESTDIR='./'
CVNREST='restart'
IFRQ_RST=0
/
&NFORCE
LINPCDF=.FALSE.
LINTERP=.FALSE.
DROFUNIT=1000
/
&NOUTPUT
COUTDIR='./'
CVARSOUT='runoff,outflw,storge'
COUTTAG='test'
IFRQ_OUT=$((hours < 24 ? hours : 24))
/
EOF

STUB_RUNOFF=$runoff "$work_dir/make/build/bin/cama.exe" cama_oasis.nml \
    > driver.log 2>&1 || { cat driver.log; exit 1; }

# One exchange per hour, discharge never negative, and river water reaching
# the ocean once the runoff has had time to flow downstream
python3 - << EOF
import numpy as np
rows = [l.split() for l in open("stub_put.log")]
sums = [float(r[3]) for r in rows]
assert len(rows) == $hours, f"expected $hours exchanges, got {len(rows)}"
assert [int(r[1]) for r in rows] == [3600 * i for i in range($hours)], rows
assert all(s >= 0 for s in sums), sums
assert sums[0] == 0.0, "no water can leave the rivers at the first exchange"
assert sums[-1] > 0.0, "no water reached the ocean: " + str(sums)
assert sums[-1] > sums[1], "outflow should increase while the runoff arrives"
print("[PASS] exchanges:", len(rows), " final discharge (kg/s): %.3e" % sums[-1])

# The runoff that CaMa-Flood sees must be the received runoff over the river
# network [kg m-2 s-1 = mm s-1] times the catchment area, in m3 s-1
nx, ny = $nx, $ny
nextx = np.fromfile("$CAMA_MAP_DIR/nextxy.bin", "<i4").reshape(2, ny, nx)[0]
area = np.fromfile("$CAMA_MAP_DIR/ctmare.bin", "<f4").reshape(ny, nx)
network = nextx != -9999
expected = float(area[network].astype(float).sum()) * $runoff / 1000.0
seen = np.fromfile("runofftest.bin", "<f4").reshape(-1, ny, nx)[0][network].sum()
assert abs(seen - expected) < 1e-3 * expected, (seen, expected)
print("[PASS] runoff into CaMa-Flood: %.4e m3/s" % seen)
EOF
grep -q "river mouths" driver.log && grep "river mouths" driver.log
echo "[PASS] CaMa-Flood OASIS driver stub test"
