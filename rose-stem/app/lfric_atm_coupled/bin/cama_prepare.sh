#!/usr/bin/env bash
##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Prepares the work directory of a coupled run in which CaMa-Flood replaces the
# JULES river model. This is run from the "cama" option of the app, before
# run_model.
#
#  1. Makes the CaMa-Flood driver available to the MOCI coupled drivers. Files
#     already present in the MOCI drivers are not replaced, so this step does
#     nothing once the driver is part of MOCI.
#  2. Replaces the river exchanges of the namcouple by those of CaMa-Flood.
#  3. Assembles a directory with the remapping weights of the coupled model
#     and those of CaMa-Flood, which MOCI links into the work directory.
#  4. Writes the CaMa-Flood dimension file for runoff on its own grid.
##############################################################################
set -euo pipefail

: "${APPS_ROOT_DIR:?}" "${CAMA_NX:?}" "${CAMA_NY:?}" "${CAMA_NLFP:?}"
: "${CAMA_WEST:?}" "${CAMA_EAST:?}" "${CAMA_NORTH:?}" "${CAMA_SOUTH:?}"
: "${CAMA_DATA:?}" "${CAMA_RMP_BASE:?}" "${CAMA_RMP_DIR:?}"
: "${CAMA_RMP_MERGED:?}"

interface_dir=$APPS_ROOT_DIR/interfaces/cama_interface

# 1. The MOCI driver. Mirror the choice of directory made by run_model.
driver_extract_dir=${DRIVER_EXTRACT_DIR:-$CYLC_WORKFLOW_SHARE_DIR/fcm_make_drivers}
if [[ $driver_extract_dir == *fcm_make* ]]; then
    drivers_dir=$driver_extract_dir/extract/drivers/Coupled_Drivers
else
    drivers_dir=$driver_extract_dir/Coupled_Drivers
fi
if [[ ! -d $drivers_dir ]]; then
    echo "[FAIL] The MOCI coupled drivers were not found in $drivers_dir" >&2
    exit 1
fi
cp -rn "$interface_dir"/moci/Coupled_Drivers/. "$drivers_dir"/

# 2. The namcouple. Everything from the river exchanges to the end of the file
# is replaced; the banner of the river section is on the 2nd of 3 lines.
snippet=namcouple_cama_rivers
n_banners=$(grep -c 'ATMOS TO RIVER' namcouple)
if [[ $n_banners != 1 ]]; then
    echo "[FAIL] Expected one river section in the namcouple, found $n_banners" >&2
    exit 1
fi
banner_line=$(grep -n 'ATMOS TO RIVER' namcouple | cut -d: -f1)
head -n $((banner_line - 2)) namcouple > namcouple.cama
sed -e "s/@CAMA_NX@/$CAMA_NX/g" -e "s/@CAMA_NY@/$CAMA_NY/g" "$snippet" \
    >> namcouple.cama
mv namcouple.cama namcouple

# 3. The remapping weights
mkdir -p "$CAMA_RMP_MERGED"
for dir in "$CAMA_RMP_BASE" "$CAMA_RMP_DIR"; do
    if [[ ! -d $dir ]]; then
        echo "[FAIL] The remapping weights directory $dir does not exist" >&2
        exit 1
    fi
    for file in "$dir"/rmp_*; do
        ln -sf "$file" "$CAMA_RMP_MERGED/$(basename "$file")"
    done
done
weights=rmp_lfricT_to_camaT_conservative1st_dstarea_nomask.nc
if [[ ! -e $CAMA_RMP_MERGED/$weights ]]; then
    echo "[FAIL] The weights $weights are not in $CAMA_RMP_DIR" >&2
    exit 1
fi

# 4. The CaMa-Flood dimension file. The runoff is supplied on the CaMa-Flood
# grid, so there is no input matrix.
cat > cama_diminfo.txt << EOF
$CAMA_NX     !! nXX
$CAMA_NY     !! nYY
$CAMA_NLFP     !! floodplain layer
$CAMA_NX     !! input nXX
$CAMA_NY     !! input nYY
1     !! input num
NONE
$CAMA_WEST     !! west  edge
$CAMA_EAST     !! east  edge
$CAMA_NORTH     !! north edge
$CAMA_SOUTH     !! south edge
EOF

mkdir -p "$CAMA_DATA"
