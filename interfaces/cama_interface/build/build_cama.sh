#!/usr/bin/env bash
##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Builds the CaMa-Flood river model as an OASIS3-MCT component (cama.exe).
#
# CaMa-Flood is not built with FCM, so this script compiles the CaMa-Flood
# source in dependency order and links it with the OASIS driver found in
# interfaces/cama_interface/source.
#
# Environment variables:
#   CAMA_SOURCE      Root of the CaMa-Flood_v4 source (required)
#   CAMA_MAKE_DIR    Install location. The executable is written to
#                    $CAMA_MAKE_DIR/build/bin/cama.exe
#                    (default: $CYLC_WORKFLOW_SHARE_DIR/$CYLC_TASK_NAME)
#   FC               Fortran compiler (default: ftn)
#   CAMA_COMPILER    cce, gnu or intel; the flags to use for FC
#                    (default: cce when FC=ftn, otherwise gnu)
#   CAMA_OMP         true to build with OpenMP (default: false)
#   OASIS_ROOT       Root of the OASIS3-MCT installation, as used by the
#                    JULES build (default: $prism_path)
#   OASIS_INC        Directories holding mod_oasis.mod
#                    (default: derived from OASIS_ROOT)
#   OASIS_LIBS       Link flags for OASIS3-MCT
#                    (default: derived from OASIS_ROOT)
#   OASIS_EXTRA_LIBS Further link flags, e.g. netCDF if FC is not a wrapper
#
# CaMa-Flood is compiled in double precision, without MPI and without netCDF;
# the coupled configuration reads and writes plain binary files only.
##############################################################################
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
interface_dir=$(cd "$script_dir/.." && pwd)

: "${CAMA_SOURCE:?CAMA_SOURCE must point at the CaMa-Flood_v4 source}"
FC=${FC:-ftn}
if [[ -z ${CAMA_COMPILER:-} ]]; then
    if [[ $FC == ftn ]]; then CAMA_COMPILER=cce; else CAMA_COMPILER=gnu; fi
fi
CAMA_OMP=${CAMA_OMP:-false}
CAMA_MAKE_DIR=${CAMA_MAKE_DIR:-${CYLC_WORKFLOW_SHARE_DIR:-.}/${CYLC_TASK_NAME:-cama}}

OASIS_ROOT=${OASIS_ROOT:-${prism_path:-}}
if [[ -z ${OASIS_INC+x} ]]; then
    : "${OASIS_ROOT:?Set OASIS_ROOT (or prism_path), or OASIS_INC and OASIS_LIBS}"
    OASIS_INC="$OASIS_ROOT/build/lib/mct $OASIS_ROOT/build/lib/psmile.MPI1 $OASIS_ROOT/build/lib/scrip"
fi
if [[ -z ${OASIS_LIBS+x} ]]; then
    : "${OASIS_ROOT:?Set OASIS_ROOT (or prism_path), or OASIS_INC and OASIS_LIBS}"
    OASIS_LIBS="-L$OASIS_ROOT/lib -lpsmile.MPI1 -lscrip -lmct -lmpeu"
fi
OASIS_EXTRA_LIBS=${OASIS_EXTRA_LIBS:-}

build_dir="$CAMA_MAKE_DIR/build"
obj_dir="$build_dir/obj"
src_dir="$build_dir/src"
bin_dir="$build_dir/bin"

# Compiler specific flags: language, module directory and OpenMP
case $CAMA_COMPILER in
    cce)
        flags="-O2 -N 255 -J $obj_dir"
        if [[ $CAMA_OMP == true ]]; then flags+=" -h omp"; else flags+=" -h noomp"; fi
        ;;
    gnu)
        flags="-O2 -cpp -ffree-line-length-none -J$obj_dir"
        if [[ $CAMA_OMP == true ]]; then flags+=" -fopenmp"; fi
        ;;
    intel)
        flags="-O2 -fpp -extend-source 400 -module $obj_dir"
        if [[ $CAMA_OMP == true ]]; then flags+=" -qopenmp"; else flags+=" -qno-openmp"; fi
        ;;
    *)
        echo "[FAIL] Unknown CAMA_COMPILER '$CAMA_COMPILER'" >&2
        exit 1
        ;;
esac
flags+=" ${CAMA_FFLAGS_EXTRA:-}"

# The CaMa-Flood modules in dependency order, as listed in src/Makefile
cama_modules="
parkind1 yos_cmf_input yos_cmf_time yos_cmf_map yos_cmf_prog yos_cmf_diag
cmf_utils_mod
cmf_calc_outflw_mod cmf_calc_pthout_mod cmf_calc_fldstg_mod
cmf_calc_stonxt_mod cmf_calc_diag_mod cmf_opt_outflw_mod
cmf_ctrl_mpi_mod cmf_ctrl_damout_mod cmf_ctrl_tracer_mod cmf_ctrl_levee_mod
cmf_ctrl_forcing_mod cmf_ctrl_boundary_mod cmf_ctrl_upsinf_mod
cmf_ctrl_output_mod cmf_ctrl_restart_mod cmf_ctrl_physics_mod
cmf_ctrl_time_mod cmf_ctrl_maps_mod cmf_ctrl_vars_mod cmf_ctrl_nmlist_mod
cmf_drv_control_mod cmf_drv_advance_mod"

if [[ ! -f $CAMA_SOURCE/src/cmf_drv_control_mod.F90 ]]; then
    echo "[FAIL] $CAMA_SOURCE does not contain the CaMa-Flood v4 source" >&2
    exit 1
fi

rm -rf "$build_dir"
mkdir -p "$obj_dir" "$bin_dir" "$src_dir"

# Compile from a private copy: compilers search the directory of a source file
# for module files, so stale build products in the shared source tree must
# not be picked up
cp "$CAMA_SOURCE"/src/*.F90 "$src_dir"
cp "$interface_dir"/source/*.F90 "$src_dir"

oasis_inc_flags=""
for dir in $OASIS_INC; do oasis_inc_flags+=" -I$dir"; done

objects=""
# The flags are lists of words, so they must not be quoted
# shellcheck disable=SC2086
for module in $cama_modules; do
    echo "[INFO] Compiling $module"
    $FC $flags -I"$obj_dir" -c "$src_dir/$module.F90" -o "$obj_dir/$module.o"
    objects+=" $obj_dir/$module.o"
done

# shellcheck disable=SC2086
for source in cama_oasis_mod cama_oasis_driver; do
    echo "[INFO] Compiling $source"
    $FC $flags -I"$obj_dir" $oasis_inc_flags \
        -c "$src_dir/$source.F90" -o "$obj_dir/$source.o"
    objects+=" $obj_dir/$source.o"
done

echo "[INFO] Linking cama.exe"
ldflags=""
case $CAMA_COMPILER in
    cce)   if [[ $CAMA_OMP == true ]]; then ldflags="-h omp"; else ldflags="-h noomp"; fi ;;
    gnu)   if [[ $CAMA_OMP == true ]]; then ldflags="-fopenmp"; fi ;;
    intel) if [[ $CAMA_OMP == true ]]; then ldflags="-qopenmp"; fi ;;
esac
# shellcheck disable=SC2086
$FC $ldflags -o "$bin_dir/cama.exe" $objects $OASIS_LIBS $OASIS_EXTRA_LIBS

echo "[INFO] Built $bin_dir/cama.exe"
