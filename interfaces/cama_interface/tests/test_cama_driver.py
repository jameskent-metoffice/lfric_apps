##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Tests for the MOCI CaMa-Flood driver (moci/Coupled_Drivers/cama_driver.py).

The driver is built on the MOCI coupled drivers library, so the tests need a
checkout of MOCI: set MOCI_COUPLED_DRIVERS to its Coupled_Drivers directory.
The tests are skipped otherwise.
"""
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

MOCI_DIR = os.environ.get("MOCI_COUPLED_DRIVERS")
CAMA_DRIVER_DIR = Path(__file__).resolve().parents[1] / "moci" / "Coupled_Drivers"

NAMELIST = """&NSIMTIME
syear=2000,
smon=1,
sday=1,
shour=0,
eyear=2001,
emon=1,
eday=1,
ehour=0,
/
&NRUNVER
lrestart=.false.,
/
&NRESTART
creststo='',
/
"""


class Environment(dict):
    """Stand-in for the common environment of the drivers."""


@unittest.skipUnless(MOCI_DIR, "MOCI_COUPLED_DRIVERS is not set")
class TestCamaDriver(unittest.TestCase):
    """The CaMa-Flood driver, run on a copy of the MOCI drivers library."""

    @classmethod
    def setUpClass(cls):
        cls.library = tempfile.mkdtemp()
        shutil.copytree(MOCI_DIR, cls.library, dirs_exist_ok=True)
        shutil.copytree(CAMA_DRIVER_DIR, cls.library, dirs_exist_ok=True)
        sys.path.insert(0, cls.library)
        # pylint: disable=import-outside-toplevel,import-error
        import cama_driver
        cls.driver = cama_driver

    @classmethod
    def tearDownClass(cls):
        sys.path.remove(cls.library)
        shutil.rmtree(cls.library)

    def setUp(self):
        self.work = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.work)
        self.previous_dir = os.getcwd()
        os.chdir(self.work)
        self.addCleanup(os.chdir, self.previous_dir)
        Path("input_cmf.nam").write_text(NAMELIST)
        self.environ = {
            "CAMA_EXEC": "/bin/true",
            "CAMA_DATA": str(Path(self.work) / "data" / "cama"),
            "ROSE_LAUNCHER_PREOPTS_CAMA": "-n 1 -ppn 1",
        }
        for name, value in self.environ.items():
            os.environ[name] = value
        self.addCleanup(self.clear_environment)
        self.common = Environment(
            CALENDAR="360day",
            TASKSTART="1979,1,1,0,0",
            TASKLENGTH="0,0,0,12,0,0",
        )

    def clear_environment(self):
        """Remove the variables that the tests set."""
        for name in list(self.environ) + ["CAMA_START"]:
            os.environ.pop(name, None)

    def run_driver(self):
        """Run the driver for a new run."""
        return self.driver.run_driver(self.common, "run_driver", {"l_namcouple": True})

    def test_new_run_sets_dates_and_starts_from_zero_storage(self):
        envar, launch, _, snd_list = self.run_driver()
        text = Path("input_cmf.nam").read_text()
        for expected in ("syear=1979,", "smon=1,", "sday=1,", "shour=0,",
                         "eyear=1979,", "emon=1,", "eday=1,", "ehour=12,",
                         "lrestart=.false.,", "creststo='',"):
            self.assertIn(expected, text)
        self.assertEqual(launch, "-n 1 -ppn 1 ./cama.exe cama_oasis.nml")
        self.assertEqual(envar["ROSE_LAUNCHER_PREOPTS_CAMA"], "'-n 1 -ppn 1'")
        self.assertTrue(os.path.islink("cama.exe"))
        self.assertTrue(Path(self.environ["CAMA_DATA"]).is_dir())
        self.assertIsNone(snd_list)

    def test_restart_file_is_used(self):
        restart = Path(self.work) / "restart1979010100.bin"
        restart.write_bytes(b"")
        os.environ["CAMA_START"] = str(restart)
        self.run_driver()
        text = Path("input_cmf.nam").read_text()
        self.assertIn("lrestart=.true.,", text)
        self.assertIn("creststo='%s'," % restart, text)

    def test_missing_restart_file_fails(self):
        os.environ["CAMA_START"] = str(Path(self.work) / "missing.bin")
        with self.assertRaises(SystemExit):
            self.run_driver()

    def test_run_crossing_a_month_end_with_a_360_day_calendar_fails(self):
        self.common["TASKSTART"] = "1979,1,20,0,0"
        self.common["TASKLENGTH"] = "0,0,10,0,0,0"
        with self.assertRaises(SystemExit):
            self.run_driver()

    def test_gregorian_run_may_cross_a_month_end(self):
        self.common["CALENDAR"] = "gregorian"
        self.common["TASKSTART"] = "1979,1,20,0,0"
        self.common["TASKLENGTH"] = "0,0,15,0,0,0"
        self.run_driver()
        text = Path("input_cmf.nam").read_text()
        self.assertIn("emon=2,", text)
        self.assertIn("eday=4,", text)

    def test_dynamic_namcouple_is_not_supported(self):
        with self.assertRaises(SystemExit):
            self.driver.run_driver(self.common, "run_driver", {"l_namcouple": False})


if __name__ == "__main__":
    unittest.main()
