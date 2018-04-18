import time
from netCDF4 import datetime
import configM2R
import model2roms
import IOstation
import clim2bry
import decimateGrid
import atmosForcing

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2018, 4, 5)
__version__ = "1.5"
__status__ = "Development"


def main():
    print("Started model2roms")
    confM2R = configM2R.Model2romsConfig()

    if confM2R.createatmosforcing or confM2R.createoceanforcing:

        if (confM2R.useesmf):
            import ESMF
            print("Starting logfile for ESMF")
           # manager = ESMF.Manager(debug=False)

        if confM2R.createoceanforcing:
          #  model2roms.convertMODEL2ROMS(confM2R)

            clim2bry.writeBry(confM2R)

        if confM2R.createatmosforcing:
            atmosForcing.createAtmosFileUV(confM2R)

    if confM2R.decimategridfile:
        decimateGrid.createGrid(confM2R.grdROMS, "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015.nc",
                                "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015v2.nc", 2)

    if confM2R.extractStations:
        print("Running in station mode and extracting pre-defined station locations")
        IOstation.getStationData(confM2R)

    print('Finished ' + time.ctime(time.time()))

main()
