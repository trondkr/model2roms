import time
from datetime import datetime
import model2roms
import IOstation
import clim2bry
import DecimateGrid
import grd
import numpy as np


__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@imr.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2014, 4, 7)
__version__ = "1.5"
__status__ = "Development"


def myhelp():
    """
    This program is run by typing: python main.py in the command window.
    """


def showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF):
    if isClimatology:
        print 'Conversions run for climatological months'
    else:
        print 'Conversions run from %s to year %s' % (start_year, end_year)
    print 'The following variables will be converted:'
    for myvar in myvars:
        print '---> %s' % myvar
    if (useESMF):
        print "All horisontal interpolations will be done using ESMF-ESMPy (module ESMF)"
    print '\nOutput grid file is: %s' % romsgridpath
    print '\nInitializing'


def main():
    print '\n--------------------------\n'
    print 'Started ' + time.ctime(time.time())

    # EDIT ===================================================================
    # Set show_progress to "False" if you do not want to see the progress
    # indicator for horizontal interpolation.
    show_progress = True
    # Set compileAll to True if you want automatic re-compilation of all the
    # fortran files necessary to run soda2roms. You need to edit compile.py for this
    compileAll = False

    # Extract time-series of data for given longitude/latitude
    extractStations = False
     # Define a set of longitude/latitude positions with names to extract into
    # station files (using extractStations)
    if (extractStations):
        stationNames = ['NorthSea', 'Iceland', 'EastandWestGreenland', 'Lofoten', 'Georges Bank']
        lonlist = [2.4301, -22.6001, -47.0801, 13.3801, -67.2001]
        latlist = [54.5601, 63.7010, 60.4201, 67.5001, 41.6423]

    # Create the bry, init, and clim files for a given grid and input data
    createForcing = True
    # Create a smaller resolution grid based on your original. Decimates every second for
    # each time run
    decimateGrid = False
    # Write ice values to file (for Arctic regions)
    writeIce = True
    # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
    useESMF = True

    # Set the input data MODEL mytype
    mytype = 'SODA'
    mytype = 'SODAMONTHLY'
    #mytype = 'GLORYS2V1'
    #mytype = 'WOAMONTHLY'
    mytype = 'NORESM'

    # Define what grid type you wnat to interpolate to:
    gridtype = "NS8KM"
    gridtype = "REGSCEN"

    # Define the paths to the input data
    if mytype == 'SODA':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODA/"
    if mytype == 'SODAMONTHLY':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODAMonthly/"
        #modelpath = "/Users/trondkr/Projects/RegScen/model2roms/testdata/"
    if mytype == 'GLORYS2V1':
        modelpath = "/Volumes/MacintoshHD2/Datasets/GLOBAL_REANALYSIS_PHYS_001_009/"
    if mytype == 'NORESM':
        modelpath = "/Users/trondkr/Projects/RegScen/NRCP45AERCN_f19_g16_CLE_01/"
     #   modelpath = "/work/users/trondk/REGSCEN/NRCP45AERCN_f19_g16_CLE_01/"
    if mytype == 'WOAMONTHLY':
        modelpath = "/Users/trondkr/Projects/is4dvar/createSSS/"

    # Define the path to the grid file
    if gridtype == "NS8KM":
        romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
    if gridtype == "REGSCEN":
        romsgridpath = "/Users/trondkr/Projects/RegScen/Grid/AA_10km_grid.nc"
     #   romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
    #    romsgridpath = "/work/users/trondk/REGSCEN/GRID/AA_10km_grid.nc"

    if mytype == 'WOAMONTHLY': isClimatology = True
    else: isClimatology = False

    start_year  = 2006
    end_year    = 2006
    start_month = 1
    end_month   = 3

    startdate = datetime(start_year, start_month, 1)
    enddate   = datetime(end_year, end_month, 1)

    # Subset the input data. The more you subset the less memory is needed for calculations
    # and the faster the process is performed. The subset is initially performed in IOsubset.py
    if gridtype == "NS8KM":
        abbreviation = "nordsjoen_8km"
        minLat = 40
        maxLat = 70
        minLon = -20
        maxLon = 40

    if gridtype == "REGSCEN":
        abbreviation = "regscen"
        minLat = -50
        maxLat = 89.5
        minLon = -179
        maxLon = 180

    subset = np.zeros(4); subset[0] = minLat; subset[1] = maxLat; subset[2] = minLon; subset[3] = maxLon

    # Name of output files for CLIM, BRY, and INIT files
    climName = abbreviation + '_clim_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    initName = abbreviation + '_init_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    bryName = abbreviation + '_bry_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    if isClimatology is True:
        climName=abbreviation + '_' + str(mytype) + '_climatology.nc'

    # Define what variables to include in the forcing files
    myvars = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
    if mytype=="NORESM":
        myvars=['temperature','salinity', 'ssh', 'uvel', 'vvel','ageice','uice','vice','aice','hice','snow_thick']
      #  myvars=['ageice','uice','vice','aice','hice','snow_thick']

    # WOA only currently contains salinity and temperature
    if isClimatology==True:
        myvars = ['temperature','salinity']

    # 5 day or 30 day average files for model input files
    if mytype == 'SODA':
        aveDays = 5.0

    if mytype in ['SODAMONTHLY', 'GLORYS2V1', 'NORESM','WOAMONTHLY']:
        aveDays = 30.0

    # NO EDIT BELOW =========================================================
    if compileAll is True:
        # Compile the Fortran 90 files to Python modules
        import compile
        compile.compileAll()

    start_day_in_start_year = np.round(((startdate - datetime(startdate.year, 1, 1)).days  + 1) / aveDays)
    end_day_in_end_year = np.round(((enddate - datetime(enddate.year, 1, 1)).days + 1) / aveDays)

    years = [(int(startdate.year) + kk) for kk in range(1 + int(enddate.year) - int(startdate.year))]

    loop = int(end_day_in_end_year) - int(start_day_in_start_year)

    if int(start_day_in_start_year) == int(end_day_in_end_year):
        IDS = [int(start_day_in_start_year) +1]
    else:
        IDS = [(i + int(start_day_in_start_year) +1) for i in range(loop + 1)]

    # FIXME: this only gives the option of running all months of the year and not subset.

    #IDS=[i  + 1 for i in range(12)]

    if isClimatology==True:
        IDS=[i+1 for i in xrange(12)]
        print "Will create climatology for months: %s"%(IDS)

    # Create the grid object for the output grid
    grdROMS = grd.grdClass(romsgridpath, "ROMS", useESMF)
    grdROMS.vars=myvars
    if (useESMF):
        # initialize MPI
        import ESMF
        manager = ESMF.Manager(logkind = ESMF.LogKind.MULTI, debug = True)

    if createForcing:

        showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF)

        model2roms.convertMODEL2ROMS(years, IDS, climName, initName, modelpath, romsgridpath, myvars, show_progress,
                                         mytype, subset, isClimatology, writeIce, useESMF)

        clim2bry.writeBry(grdROMS, start_year, bryName, climName, writeIce, mytype)

    if decimateGrid:
        DecimateGrid.createGrid(grdROMS, '/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc', 2)

    if extractStations:
        print "Running in station mode and extracting pre-defined station locations"
        IOstation.getStationData(years, IDS, modelpath, latlist, lonlist, stationNames)

    print 'Finished ' + time.ctime(time.time())


if __name__ == "__main__":
    main()
