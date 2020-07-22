import unittest
from datetime import datetime

import numpy as np
import xarray as xr

import model2roms
import grd
import configM2R

class TestModel2roms(unittest.TestCase):

    def setUp(self):
        self.configM2R = configM2R.Model2romsConfig()
        self.assertIsNotNone(self.configM2R)

class TestConfigM2R(TestModel2roms):

    def test_output_names(self):
        self.assertTrue(self.configM2R.clim_name)
        self.assertTrue(self.configM2R.init_name)
        self.assertTrue(self.configM2R.bry_name)

    def test_define_ocean_forcing_data_path_returns_correct(self):
        self.configM2R.ocean_indata_type = "SODA3"
        self.assertTrue(self.configM2R.define_ocean_forcing_data_path())
        self.configM2R.ocean_indata_type = "GLORYS"
        self.assertTrue(self.configM2R.define_ocean_forcing_data_path())
        self.configM2R.ocean_indata_type = "SODA3_5DAY"
        self.assertTrue(self.configM2R.define_ocean_forcing_data_path())

    def test_define_ocean_forcing_data_path_returns_error_when_missing_dataset(self):
        self.configM2R.ocean_indata_type = "not_valid"
        self.assertEqual(self.configM2R.define_ocean_forcing_data_path(), KeyError)

if __name__ == "__main__":
    unittest.main()
