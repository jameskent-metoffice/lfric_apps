# CaMa-Flood interface

Couples the [CaMa-Flood](https://github.com/global-hydrodynamics/CaMa-Flood_v4)
river routing model to the LFRic coupled model, in place of the JULES river
model (`river.exe`). CaMa-Flood runs as an OASIS3-MCT component next to LFRic and
NEMO.

## How it works

```
LFRic ──surf_roff, sub_surf_roff──▶ CaMa-Flood ──outflow_per_river──▶ NEMO
        [kg m-2 s-1, CaMa grid]                   [kg s-1, one per river]
```

The exchanges, field names, grids and coupling period are the same as those of
`river.exe`, so LFRic (`interfaces/coupled_interface`) and NEMO do not change.
Only the river section of the namcouple does.

* LFRic sends its runoff, remapped by OASIS onto the CaMa-Flood map grid
  (`camaT`). CaMa-Flood reads it on that grid with no input matrix of its own
  (`LINTERP=.FALSE.`), in units of mm s-1 (`DROFUNIT=1000`).
* At each coupling time the driver sends the river discharge first and then
  receives the runoff, and then advances CaMa-Flood by one coupling period
  (the CaMa-Flood time step is the coupling period, 3600 s). This is the
  sequence of `river.exe` and the `SEQ` values of the namcouple depend on it.
* The discharge of every CaMa-Flood river mouth that drains to the ocean is
  added to the NEMO river number found there in a map of river numbers on the
  CaMa-Flood grid, and converted to kg s-1 (1000 kg m-3).
* Inland basins (CaMa-Flood `nextx = -10`) are not sent to the ocean, as with the
  JULES river model; the water that CaMa-Flood sends out of them is lost.

## Contents

| Path | Purpose |
| ---- | ------- |
| `source/cama_oasis_mod.F90` | OASIS coupling layer |
| `source/cama_oasis_driver.F90` | Program that runs CaMa-Flood and the exchanges |
| `build/build_cama.sh` | Builds `cama.exe` from the CaMa-Flood source and the above |
| `tools/make_cama_ancils.py` | Creates the river number map, `CDIMINFO` file and SCRIP grid |
| `moci/Coupled_Drivers/` | The MOCI driver for the `cama` component (see below) |
| `tests/` | Tests, none of which need OASIS or the Met Office data |

In `rose-stem/`:

| Path | Purpose |
| ---- | ------- |
| `app/fcm_make_cama/` | Builds CaMa-Flood (runs `build_cama.sh`) |
| `app/lfric_atm_coupled/opt/rose-app-cama.conf` | The `cama` option of the coupled app: CaMa-Flood namelists, `models`, command |
| `app/lfric_atm_coupled/bin/cama_prepare.sh` | Prepares the work directory before `run_model` |
| `app/lfric_atm_coupled/file/namcouple_cama_rivers` | The river section of the namcouple for CaMa-Flood |
| `templates/` | The `cama` task setting and its defaults, environment and build task |

CaMa-Flood itself is a dependency in `dependencies.yaml` (`CaMa-Flood_v4`, v4.31).

## Running the test

A task sets `"cama": true` in its definition and adds the `cama` option to its
`opt_confs`. The test `nwp_gc6_cama-C48_MG` does this; it is the same as
`nwp_gc6-C48_MG` otherwise. Run it with the group `lfric_atm_coupled_cama_ex1a`.
It has no known good output yet: after the first successful run set
`kgo_checks` to `["checksum"]` and add the checksum.

Tasks without `"cama": true` are unchanged and still build and run `river.exe`.

## Data that is needed

None of these exist yet. The default locations are in
`templates/default_task_definitions_coupled.cylc` and can be set per task
(`cama_map_dir`, `cama_river_number_file`, `cama_rmp_dir`, `cama_start`).

1. **CaMa-Flood river map**, for the global 15 minute map: `nextxy.bin`,
   `ctmare.bin`, `elevtn.bin`, `nxtdst.bin`, `rivlen.bin`, `fldhgt.bin`,
   `rivwth_gwdlr.bin`, `rivhgt.bin`, `rivman.bin` and `bifprm.txt`. They are part
   of the CaMa-Flood distribution, which requires registration and has its own
   licence terms. The grid is given by `cama_nx`, `cama_ny`, `cama_west` etc.
2. **River number map** `cama_river_number_file`, the NEMO river number of each
   river mouth on the CaMa-Flood grid:
   ```
   make_cama_ancils.py river-numbers --map-dir <map> --river-number-file \
       <river_number_trip.nc> --output cama_river_number.bin
   ```
   `river_number_trip.nc` is the JULES river number ancillary that NEMO's river
   numbers (`CPL_RIVER_COUNT`, 1681) come from. Check the report for mouths that
   are not mapped, which are not sent to the ocean.
3. **Remapping weights** from LFRic to the CaMa-Flood grid,
   `rmp_lfricT_to_camaT_conservative1st_dstarea_nomask.nc`, in `cama_rmp_dir`.
   `make_cama_ancils.py scrip-grid` writes a description of the CaMa-Flood grid
   for the tool that made the LFRic to river weights (`generate_weights`). The
   rows of the grid run from north to south and the destination index is
   `x + (y - 1) * nx`. `cama_prepare.sh` links these with the weights of the
   coupled model, so `cama_rmp_dir` only needs to hold the new file.
4. **A restart file** (optional), `cama_start`. Without one CaMa-Flood starts
   from zero storage (which is not empty: CaMa-Flood starts with some water in
   its rivers). Segments of a continuation run restart from the latest
   `restart*.bin` of the previous segment.

## Building

`build_cama.sh` compiles CaMa-Flood (double precision, no MPI, no netCDF: all
files are plain binary) and links it with OASIS3-MCT, found from `OASIS_ROOT` or
`prism_path` as for the JULES build. See the header of the script for the
options, in particular `FC`/`CAMA_COMPILER` (CCE, GNU and Intel are supported).
CaMa-Flood runs on a single rank (`ROSE_LAUNCHER_PREOPTS_CAMA = -n 1 -ppn 1`).

## The MOCI driver

The coupled model is launched by the MOCI `Coupled_Drivers`, which need a driver
for every entry of `models`. MOCI at the pinned version has `rivers_driver.py`
but none for CaMa-Flood, so `moci/Coupled_Drivers/` holds `cama_driver.py` and
`dr_env_lib/cama_def.py`. `cama_prepare.sh` copies them next to the MOCI drivers
without replacing existing files, so this stops being used when the driver is
added to MOCI, which is where it belongs; then these files can be removed.

The MCT driver of MOCI does not know the `cama` component, so `cama` is not in
`COUPLING_COMPONENTS` (it only uses that to clean up debug files). The static
namcouple is what couples it.

The driver sets the run dates and the restart file in the CaMa-Flood namelist.
CaMa-Flood uses the Gregorian calendar, so with the 360 day calendar of the
coupled tests it can only run within one month of 28 days or less.

## Tests

* `tests/test_make_cama_ancils.py`: `python3 -m unittest test_make_cama_ancils`
* `tests/test_cama_driver.py`: needs a MOCI checkout, for example
  `MOCI_COUPLED_DRIVERS=<moci>/Coupled_Drivers python3 -m unittest test_cama_driver`
* `tests/run_stub_test.sh`: builds `cama.exe` against stand-in OASIS and MPI
  modules (`tests/stubs`) and runs it on a CaMa-Flood map, for example
  `map/mekong_15min` of CaMa-Flood. It checks the exchanges, the mapping of river
  mouths, and that the runoff in CaMa-Flood is the runoff received.
  Needs `CAMA_SOURCE` and `CAMA_MAP_DIR`; see the script for the details.

## Not tested

The code has been run with CaMa-Flood on the Mekong map and a stand-in for OASIS,
and the suite has been validated and rendered, but not with the real OASIS3-MCT,
the CCE compiler, a global map or a coupled run. The LFRic to CaMa-Flood weights
do not exist yet.
