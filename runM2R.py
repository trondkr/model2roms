import time
from datetime import datetime
import logging
import IOstation
import clim2bry
import configM2R
import decimateGrid
import model2roms

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2021, 7, 27)
__version__ = "1.6"
__status__ = "Development"

"""
    Main method for running model2roms
    Start: python runM2R.py 
"""

def run():
    logging.basicConfig(level=logging.INFO)
    logging.info("[M2R_run] Initialized logging")

    logging.info("[M2R_run] Started model2roms")
    confM2R = configM2R.Model2romsConfig()
    confM2R.create_grd_objects()

    if confM2R.create_atmos_forcing or confM2R.create_ocean_forcing:

        if confM2R.create_ocean_forcing:
            model2roms.convert_MODEL2ROMS(confM2R)

            clim2bry.writebry(confM2R)

      #  if confM2R.createAtmosForcing:
      #      atmosForcing.createAtmosFileUV(confM2R)

    if confM2R.decimate_gridfile:
        decimateGrid.createGrid(confM2R.grdROMS, "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015.nc",
                                "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015v2.nc", 2)

    if confM2R.extract_stations:
        print("Running in station mode and extracting pre-defined station locations")
        IOstation.getStationData(confM2R)

    print('Finished ' + time.ctime(time.time()))

run()
