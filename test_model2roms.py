import unittest
from datetime import datetime

import numpy as np
import xarray as xr

import model2roms
import grd
import configM2R

# Unittest for ``CMIP6_GLORYS12V1` setup

class TestModel2roms(unittest.TestCase):
    def setUp(self):
        self.configM2R=confM2R = configM2R.Model2romsConfig()
        self.assertIsNone(self.configM2R)

if __name__ == "__main__":
    unittest.main()
