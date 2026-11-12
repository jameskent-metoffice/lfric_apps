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
    cama_def.py

DESCRIPTION
    Definition of the environment variables required for a CaMa-Flood river
    model run
'''
CAMA_ENVIRONMENT_VARS_INITIAL = {
    'CAMA_EXEC': {},
    'CAMA_LINK': {'default_val': 'cama.exe'},
    # CaMa-Flood namelist, with the run dates and restart file set by the
    # driver
    'CAMA_NLIST': {'default_val': 'input_cmf.nam'},
    # Namelist of the OASIS coupling of CaMa-Flood
    'CAMA_COUPLE_NLIST': {'default_val': 'cama_oasis.nml'},
    # Directory of the CaMa-Flood restart and output files
    'CAMA_DATA': {},
    # CaMa-Flood restart file, or unset to start from zero storage
    'CAMA_START': {'default_val': 'unset'},
    'ROSE_LAUNCHER_PREOPTS_CAMA': {},
}
