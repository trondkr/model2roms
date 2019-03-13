import time, calendar
from netCDF4 import Dataset, date2num, num2date
import model2roms
import IOstation
import clim2bry
import decimateGrid
import grd
import numpy as np
import atmosForcing
import sys
from datetime import datetime, timedelta

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2018, 4, 5)
__version__ = "1.5"
__status__ = "Development"


class Model2romsConfig(object):

    def definesubsetforindata(self):
        # Subset the input data. The more you subset the less memory is needed for calculations
        # and the faster the process is performed. The subset is initially performed in IOsubset.py
        subset = np.zeros(4)

        if self.outgrid == "NS8KM":
            subset[1] = 40
            subset[1] = 70
            subset[2] = -30
            subset[3] = 40

        if self.outgrid == "A20":
            subset[0] = 30
            subset[1] = 90
            subset[2] = -179
            subset[3] = 360

        return subset

    # Define abbreviation for the run: sued to name output files etc.
    def defineabbreviation(self):
        return {"A20": "a20",
                "ROHO800": "roho800"}[self.outgrid]

    def showinfo(self):
        if self.isclimatology:
            print('\n=>Conversions run for climatological months')
        else:
            print('\n=>Conversions run from %s to year %s' % (self.start_year, self.end_year))
        print('==> The following variables will be converted:')
        for myvar in self.globalvarnames:
            print('===> %s' % myvar)
        if self.useesmf:
            print("=>All horisontal interpolations will be done using ESMF-ESMPy (module ESMF)")
        print("=>Output files are written in format: %s" % self.myformat)
        print('\n=>Output grid file is: %s' % self.romsgridpath)

    def formatdatesforoutputnames(self):
        # Format the date for use in output filenames
        startmonth = ("0%s" % self.start_month if self.start_month < 10 else "%s" % self.start_month)
        startday = ("0%s" % self.start_day if self.start_day < 10 else "%s" % self.start_day)
        endmonth = ("0%s" % self.end_month if self.end_month < 10 else "%s" % self.end_month)
        endday = ("0%s" % self.end_day if self.end_day < 10 else "%s" % self.end_day)

        modelperiod = str(self.start_year) + str(startmonth) + str(startday) + '_to_' + str(self.end_year) + str(
            endmonth) + str(
            endday)

        return modelperiod

    def defineoutputfilenames(self):
        # Get string representation of start and end dates
        modelperiod = self.formatdatesforoutputnames()

        # Name of output files for CLIM, BRY, and INIT files
        climname = self.abbreviation + '_clim_' + str(self.indatatype) + '_' + str(modelperiod) + '.nc'
        initname = self.abbreviation + '_init_' + str(self.indatatype) + '_' + str(modelperiod) + '.nc'
        bryname = self.abbreviation + '_bry_' + str(self.indatatype) + '_' + str(modelperiod) + '.nc'

        return climname, initname, bryname

    # Define the global variables to be used for each type of input data. Not all input datasets contains information on
    # e.g. sea ice so those variables can not be included. the SODA3si for example does not contain ssh.
    # OPTIONS: ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'ageice', 'uice', 'vice', 'aice', 'hice']

    def defineglobalvarnames(self):

        return {'SODA': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel'],
                'SODA3': ['temperature', 'salinity', 'uvel', 'vvel'],
                'GLORYS': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'uice', 'vice', 'aice', 'hice'],
                'WOAMONTHLY': ['temperature', 'salinity'],
                'NORESM': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'ageice', 'uice', 'vice', 'aice', 'hice','hs',
                'O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']}[
            self.indatatype]

    # Define the corresponding name of the variables in the input dataset files. This list needs to correspond
    # exactly with the list given in the function defineglobalvarnames:

    def defineinputdatavarnames(self):

        return {'SODA': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel'],
                'SODA3': ['temp', 'salt', 'ssh', 'u', 'v'],
                'GLORYS': ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty', 'iicevelu', 'iicevelv',
                           'ileadfra', 'iicethic'],
                'WOAMONTHLY': ['temperature', 'salinity'],
                'NORESM': ['templvl', 'salnlvl', 'sealv', 'uvellvl', 'vvellvl', 'iage', 'uvel', 'vvel', 'aice', 'hi',
                           'hs','dissic','talk','po4','no3','si','o2']}[self.indatatype]

    def defineromsgridpath(self):
        return {'A20': '/cluster/projects/nn9412k/A20/Grid/A20niva_grd_v1.nc',
                'ROHO800': '/global/work/andrest/ROHO800/River_and_Grid_fix/ROHO800_grid_fix.nc'}[self.outgrid]

    def defineforcingdatapath(self):
        return {'SODA3': "/global/home/trondk/Projects/SODA3/SODA3/dsrs.atmos.umd.edu/DATA/soda3.3.1/REGRIDED/",
                'NORESM': "/cluster/projects/nn9412k/A20/FORCING/RCP85_ocean/",
                'WOAMONTHLY': "/Users/trondkr/Projects/is4dvar/createSSS/"}[self.indatatype]

    def __init__(self):
        print('\n--------------------------\n')
        print('Started ' + time.ctime(time.time()))

        # EDIT ===================================================================
        # Set showprogress to "False" if you do not want to see the progress
        # indicator for horizontal interpolation.
        self.showprogress = True
        # Set compileAll to True if you want automatic re-compilation of all the
        # fortran files necessary to run model2roms. Options are "gfortran" or "ifort". Edit
        # compile.py to add other Fortran compilers.
        self.compileall = False
        # Extract time-series of data for given longitude/latitude
        self.extractstations = False
        # Define a set of longitude/latitude positions with names to extract into
        # station files (using extractStations)
        if self.extractstations:
            #  stationNames = ['NorthSea', 'Iceland', 'EastandWestGreenland', 'Lofoten', 'Georges Bank']
            #  lonlist = [2.4301, -22.6001, -47.0801, 13.3801, -67.2001]
            #  latlist = [54.5601, 63.7010, 60.4201, 67.5001, 41.6423]

            self.stationnames = ["Ytre Utsira", "Indre Utsira", "Lista"]
            self.latlist = [59.316667, 59.316667, 58.016667]
            self.lonlist = [4.800000, 4.983333, 6.533333]
            self.numberofpoints = 4 # Number of points around lat/lon to extract and average as output

        # Create the bry, init, and clim files for a given grid and input data
        self.createoceanforcing = True
        # Create atmospheric forcing for the given grid
        self.createatmosforcing = False  # currently in beta stages and unavailable
        # Create a smaller resolution grid based on your original. Decimates every second for
        # each time run
        self.decimategridfile = False
        # Write ice values to file (for Arctic regions)
        self.writeice = True
        # Write biogeochemistry values to file
        self.writebcg = True
        # ROMS sometimes requires input of ice and ssh, but if you dont have these write zero files to file
        self.set2DvarsToZero=True
        # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
        self.useesmf = True
        # Apply filter to smooth the 2D fields after interpolation (time consuming but enhances results)
        self.usefilter = True
        # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
        # Using NETCDF4 automatically turns on compression of files (ZLIB)
        self.myformat = 'NETCDF4'
        self.myzlib = True
        # Frequency of the input data: usually monthly
        self.timefrequencyofinputdata = "month"  # , "month", "hour"

        # IN GRIDTYPES ------------------------------------------------------------------------------
        #  Define what grid type you wnat to interpolate from (input MODEL data)
        # Options:
        # 1. SODA, 2. SODAMONTHLY, 3.WOAMONTHLY, 4. NORESM, 4. GLORYS, 5. SODA3
        self.indatatype = 'NORESM'

        # Define contact info for final NetCDF files
        self.authorname = "Trond Kristiansen"
        self.authoremail = "trond.kristiansen (at) niva.no"

        # Define what grid type you wnat to interpolate from: Can be Z for SIGMA for ROMS
        # vertical coordinate system or ZLEVEL. also define the name of the dimensions in the input files.
        # Options:
        # 1. SIGMA (not prpoerly implemented yet), 2. ZLEVEL
        self.ingridtype = "ZLEVEL"

        # Define the names of the geographical variables in the input files
        self.grdtype = 'regular'
        self.lonname = "plon"
        self.latname = "plat"
        self.depthname = "pdepth"
        self.lonname_u = "ulon"
        self.latname_u = "ulat"
        self.lonname_v = "vlon"
        self.latname_v = "vlat"
        self.timename = "time"
        self.realm = "ocean"
        self.fillvaluein = -1.e20

        # OUT GRIDTYPES ------------------------------------------------------------------------------
        # Define what grid type you wnat to interpolate to
        # Options: This is just the name of your grid used to identify your selection later
        self.outgrid = "A20" #"ROHO800"
        self.outgridtype = "ROMS"

        # Subset input data. If you have global data you may want to seubset these to speed up reading. Make
        # sure that your input data are cartesian (0-360 or -180:180, -90:90)
        self.subsetindata = False
        if self.subsetindata:
            self.subset = self.definesubsetforindata()

        # Define nmber of output depth levels
        self.nlevels = 40
        # Define the grid stretching properties (leave default if uncertain what to pick)
        self.vstretching = 4
        self.vtransform = 2
        self.theta_s = 7.0
        self.theta_b = 0.1
        self.tcline = 250.0
        self.hc = 250

        # PATH TO FORCINGDATA --------------------------------------------------------------------
        # Define the path to the input data
        self.modelpath = self.defineforcingdatapath()

        # PATH TO GRID -----------------------------------------------------------------------------
        # Define the path to the grid file
        self.romsgridpath = self.defineromsgridpath()

        # Climatology is only monthly and model2roms needs to know this
        self.isclimatology = True if self.indatatype == 'WOAMONTHLY' else False

        # DATE AND TIME DETAILS ---------------------------------------------------------
        # Define the period to create forcing for
        self.start_year = 2006
        self.end_year = 2007
        self.start_month = 1
        self.end_month = 3
        self.start_day = 15
        self.end_day = 15

        if int(calendar.monthrange(self.start_year, self.start_month)[1]) < self.start_day:
            self.start_day = int(calendar.monthrange(self.start_year, self.start_month)[1])

        if int(calendar.monthrange(self.end_year, self.end_month)[1]) < self.end_day:
            self.end_day = int(calendar.monthrange(self.end_year, self.end_month)[1])

        self.startdate = datetime(self.start_year, self.start_month, self.start_day)
        self.enddate = datetime(self.end_year, self.end_month, self.end_day)
        self.years = [self.start_year + year for year in range(self.end_year + 1 - self.start_year)]

        # DEFINE VARIABLE NAMES ---------------------------------------------------------
        # Define what and name of variables to include in the forcing files
        # -> myvars is the name model2roms uses to identify variables
        # -> varNames is the name of the variable found in the NetCDF input files

        self.globalvarnames = self.defineglobalvarnames()
        self.inputdatavarnames = self.defineinputdatavarnames()

        # NO EDIT BELOW ====================================================================================================
        if self.compileall is True:
            import compile;
            compile.compileallgfortran()

        if self.createatmosforcing or self.createoceanforcing:
            self.abbreviation = self.defineabbreviation()

            self.climname, self.initname, self.bryname = self.defineoutputfilenames()

            if self.isclimatology is True:
                self.climname = self.abbreviation + '_' + str(self.indatatype) + '_climatology.nc'

            self.showinfo()

            if self.useesmf:
                try:
                    import ESMF
                except ImportError:
                    raise ImportError("Unable to import ESMF")
                print("Starting logfile for ESMF")
                ESMF.Manager(debug=True)

            # Create the grid object for the output grid
            self.grdROMS = grd.Grd("ROMS", self)
            self.grdROMS.nlevels =self.nlevels
            self.grdROMS.vstretching = self.vstretching
            self.grdROMS.vtransform = self.vtransform
            self.grdROMS.theta_s = self.theta_s
            self.grdROMS.theta_b = self.theta_b
            self.grdROMS.tcline = self.tcline
            self.grdROMS.hc = self.hc
            self.grdROMS.lonname = 'lon_rho'
            self.grdROMS.latname = 'lat_rho'

            self.grdROMS.opennetcdf(self.romsgridpath)
            self.grdROMS.createobject(self)
            self.grdROMS.getdims()

            # Create the grid object for the input grid
            self.grdMODEL = grd.Grd("FORCINGDATA", self)
            self.grdMODEL.grdType = self.grdtype
            self.grdMODEL.lonName = self.lonname
            self.grdMODEL.latName = self.latname
            self.grdMODEL.depthName = self.depthname
            self.grdMODEL.fillval = self.fillvaluein
