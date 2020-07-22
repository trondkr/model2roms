import unittest
from datetime import datetime

import numpy as np
import xarray as xr

import model2roms

# Unittest for ``CMIP6_GLORYS12V1` setup

class TestModel2roms(unittest.TestCase):
    def setUp(self):
        self.assertTrue(model2roms)


class TestMethods(TestModel2roms):
    def test_init(self):
        self.assertTrue(model2roms)

if __name__ == "__main__":
    unittest.main()
