#!/usr/bin/env python3
'''
*****************************COPYRIGHT******************************
 (C) Crown copyright 2026 Met Office. All rights reserved.

 Use, duplication or disclosure of this code is subject to the restrictions
 as set forth in the licence. If no licence has been raised with this copy
 of the code, the use, duplication or disclosure of it is strictly
 prohibited. Permission to do so must first be obtained in writing from the
 Met Office Information Asset Owner at the following address:

 Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB, United Kingdom
*****************************COPYRIGHT******************************
NAME
    cama_driver.py

DESCRIPTION
    Driver for the CaMa-Flood river model, called from link_drivers. This
    takes the place of rivers_driver.py when CaMa-Flood is coupled to LFRic
    and NEMO instead of the JULES river model. The coupling is described by
    a static namcouple, so nothing is generated at run time.
'''
import datetime
import os
import pathlib
import sys
import common
import error
import dr_env_lib.cama_def
import dr_env_lib.env_lib

# The longest run for which the calendar of the suite and the Gregorian
# calendar of CaMa-Flood can not disagree about the end date: the run must
# stay within a month of 28 days or more
MAX_NON_GREGORIAN_DAYS = 28


def _end_of_run(common_envar):
    '''
    Find the start and end date of the run. CaMa-Flood uses the Gregorian
    calendar, so other calendars can only be used for runs that do not
    cross the end of a month.
    '''
    calendar = common_envar['CALENDAR']
    if calendar not in ('360day', '365day', 'gregorian'):
        sys.stderr.write('[FAIL] Calendar type %s not recognised\n' %
                         calendar)
        sys.exit(error.INVALID_EVAR_ERROR)

    start = [int(i) for i in common_envar['TASKSTART'].split(',')]
    start_date = datetime.datetime(start[0], start[1], start[2], start[3],
                                   start[4])
    run_seconds = common.setup_runtime(common_envar)
    end_date = start_date + datetime.timedelta(seconds=run_seconds)

    if calendar != 'gregorian':
        if (start_date.day - 1) * 86400 + run_seconds >= \
                MAX_NON_GREGORIAN_DAYS * 86400:
            sys.stderr.write('[FAIL] CaMa-Flood uses the Gregorian calendar. '
                             'With the %s calendar it can only run within a '
                             'month of %i days.\n' %
                             (calendar, MAX_NON_GREGORIAN_DAYS))
            sys.exit(error.INVALID_EVAR_ERROR)
    if end_date.minute or end_date.second or start_date.minute:
        sys.stderr.write('[FAIL] CaMa-Flood can only start and end a run on '
                         'the hour\n')
        sys.exit(error.INVALID_EVAR_ERROR)

    return start_date, end_date


def _update_cama_nl(cama_envar, start_date, end_date):
    '''
    Check that the CaMa-Flood namelist exists, set the start and end dates
    and the restart file, and create the data directory
    '''
    nlist = cama_envar['CAMA_NLIST']
    if not os.path.isfile(nlist):
        sys.stderr.write('[FAIL] Can not find the CaMa-Flood namelist file '
                         '%s\n' % nlist)
        sys.exit(error.MISSING_DRIVER_FILE_ERROR)

    mod_nlist = common.ModNamelist(nlist)
    for prefix, date in (('s', start_date), ('e', end_date)):
        mod_nlist.var_val(prefix + 'year', date.year)
        mod_nlist.var_val(prefix + 'mon', date.month)
        mod_nlist.var_val(prefix + 'day', date.day)
        mod_nlist.var_val(prefix + 'hour', date.hour)

    # Without a restart file CaMa-Flood starts from zero storage
    restart = cama_envar['CAMA_START']
    if restart in ('', 'unset'):
        mod_nlist.var_val('lrestart', '.false.')
        mod_nlist.var_val('creststo', '')
    else:
        if not os.path.isfile(restart):
            sys.stderr.write('[FAIL] The CaMa-Flood restart file %s does not '
                             'exist\n' % restart)
            sys.exit(error.RESTART_FILE_ERROR)
        sys.stdout.write('[INFO] CaMa-Flood restarts from %s\n' % restart)
        mod_nlist.var_val('lrestart', '.true.')
        mod_nlist.var_val('creststo', restart)
    mod_nlist.replace()

    # CaMa-Flood writes its restart and output files here
    pathlib.Path(cama_envar['CAMA_DATA']).mkdir(parents=True, exist_ok=True)


def _setup_executable(common_envar):
    '''
    Setup the environment and any files required by the executable
    '''
    cama_envar = dr_env_lib.env_lib.LoadEnvar()
    cama_envar = dr_env_lib.env_lib.load_envar_from_definition(
        cama_envar, dr_env_lib.cama_def.CAMA_ENVIRONMENT_VARS_INITIAL)

    # Link the executable
    common.remove_file(cama_envar['CAMA_LINK'])
    os.symlink(cama_envar['CAMA_EXEC'], cama_envar['CAMA_LINK'])

    start_date, end_date = _end_of_run(common_envar)
    _update_cama_nl(cama_envar, start_date, end_date)

    return cama_envar


def _set_launcher_command(cama_envar):
    '''
    Setup the launcher command for the executable. The coupling configuration
    namelist is passed as the argument.
    '''
    launch_cmd = '%s ./%s %s' % (cama_envar['ROSE_LAUNCHER_PREOPTS_CAMA'],
                                 cama_envar['CAMA_LINK'],
                                 cama_envar['CAMA_COUPLE_NLIST'])

    # Put in quotes to allow this environment variable to be exported as it
    # contains (or can contain) spaces
    cama_envar['ROSE_LAUNCHER_PREOPTS_CAMA'] = "'%s'" % \
        cama_envar['ROSE_LAUNCHER_PREOPTS_CAMA']
    return launch_cmd


def _finalize_executable():
    '''
    Finalize the CaMa-Flood run. The restart file is written by CaMa-Flood
    to CAMA_DATA, from where the next cycle finds it.
    '''
    sys.stdout.write('[INFO] finalizing CaMa-Flood\n')


def run_driver(common_env, mode, run_info):
    '''
    Run the driver, and return an instance of common.LoadEnvar and a string
    containing the launcher command for the CaMa-Flood model
    '''
    model_snd_list = None
    if mode == 'run_driver':
        if not run_info['l_namcouple']:
            sys.stderr.write('[FAIL] The CaMa-Flood driver needs a static '
                             'namcouple file\n')
            sys.exit(error.INVALID_DRIVER_ARG_ERROR)
        exe_envar = _setup_executable(common_env)
        launch_cmd = _set_launcher_command(exe_envar)
    elif mode == 'finalize':
        _finalize_executable()
        exe_envar = None
        launch_cmd = None
    return exe_envar, launch_cmd, run_info, model_snd_list
