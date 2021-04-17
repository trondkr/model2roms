import calendar
import logging
import os
import time
from datetime import datetime

import numpy as np

import grd

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2021, 3, 23)
__version__ = "1.6"
__status__ = "Development"


class Model2romsConfig(object):

    def define_subset_for_indata(self):
        # Subset the input data. The more you subset the less memory is needed for calculations
        # and the faster the process is performed. The subset is initially performed in IOsubset.py
        subset = np.zeros(4)

        if self.outgrid_name == "NS8KM":
            return subset[40, 70, -30, 40]

        elif self.outgrid_name == "A20":
            return subset[30, 90, -179, 360]
        else:
            raise Exception("Unable to subset {}".format(self.outgrid_name))

    def showinfo(self):
        if self.isclimatology:
            logging.info('[M2R_configM2R]\n=> Conversions run for climatological months')
        else:
            logging.info('[M2R_configM2R]\n=> Conversions run from year/month: %s/%s to %s/%s' % (
                self.start_year, self.start_month, self.end_year, self.end_month))
        logging.info('[M2R_configM2R]==> The following variables will be interpolated: {}'.format(self.global_varnames))

        if self.use_esmf:
            logging.info('[M2R_configM2R]=>All horisontal interpolations will be done using ESMF')
        logging.info('[M2R_configM2R] => Output files are written in format: {}'.format(self.output_format))
        logging.info('[M2R_configM2R] => Output grid file is: {}'.format(self.roms_grid_path))

    def format_dates_for_outputnames(self) -> str:
        return "{}{}{}_to_{}{}{}".format(str(self.start_year),
                                         str(self.start_month).zfill(2),
                                         str(self.start_day).zfill(2),
                                         str(self.end_year).zfill(2),
                                         str(self.end_month).zfill(2),
                                         str(self.end_day).zfill(2))

    def define_output_filenames(self) -> (str, str, str):
        # Get string representation of start and end dates
        modelperiod = self.format_dates_for_outputnames()

        # Name of output files for CLIM, BRY, and INIT files
        climname = self.abbreviation + '_clim_' + str(self.ocean_indata_type) + '_' + str(modelperiod) + '.nc'
        initname = self.abbreviation + '_init_' + str(self.ocean_indata_type) + '_' + str(modelperiod) + '.nc'
        bryname = self.abbreviation + '_bry_' + str(self.ocean_indata_type) + '_' + str(modelperiod) + '.nc'

        return climname, initname, bryname

    # Define the global variables to be used for each type of input data. Not all input datasets contains information on
    # e.g. sea ice so those variables can not be included. the SODA3si for example does not contain ssh.
    # OPTIONS: ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'ageice', 'uice', 'vice', 'aice', 'hice']

    def define_global_varnames(self):

        return {'SODA': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel'],
                'SODA3': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel'],
                'SODA3_5DAY': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel'],
                'GLORYS': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'uice', 'vice', 'aice', 'hice'],
                'WOAMONTHLY': ['temperature', 'salinity'],
                'NORESM': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'ageice', 'uice', 'vice', 'aice', 'hice',
                           'hs',
                           'O3_c', 'O3_TA', 'N1_p', 'N3_n', 'N5_s', 'O2_o']}[
            self.ocean_indata_type]

    # Define the corresponding name of the variables in the input dataset files. This list needs to correspond
    # exactly with the list given in the function define_global_varnames:

    def define_input_data_varnames(self):
        return {'SODA3': ['temp', 'salt', 'ssh', 'u', 'v'],
                'SODA3_5DAY': ['temp', 'salt', 'ssh', 'u', 'v'],
                'GLORYS': ['thetao', 'so', 'zos', 'uo', 'vo', 'usi', 'vsi',
                           'siconc', 'sithick'],
                'NORESM': ['templvl', 'salnlvl', 'sealv', 'uvellvl', 'vvellvl', 'iage', 'uvel', 'vvel', 'aice', 'hi',
                           'hs', 'dissic', 'talk', 'po4', 'no3', 'si', 'o2']}[self.ocean_indata_type]

    # Define the path to where the  ROMS grid can be found
    def define_roms_grid_path(self):
        try:
            return {'Antarctic': '../oceanography/model2roms_grids/ANT_01.nc',
                    'A20': '/Users/trondkr/Dropbox/NIVA/A20/Grid/A20niva_grd_v1.nc',
                    # '/cluster/projects/nn9412k/A20/Grid/A20niva_grd_v1.nc',
                    #     'ROHO800': '/cluster/projects/nn9490k/ROHO800/Grid/ROHO800_grid_fix3.nc'}[self.outgrid_name]
                    'ROHO800': '/Users/trondkr/Dropbox/NIVA/ROHO800/Grid/ROHO800_grid_fix3.nc'}[self.outgrid_name]
        except KeyError:
            return KeyError

    # Define the abbreviation for the run, which is used to name output files etc.
    def define_abbreviation(self):
        return {"A20": "a20",
                "Antarctic": "Antarctic",
                "ROHO800": "roho800"}[self.outgrid_name]

    def define_ocean_forcing_data_path(self):
        try:
            return {'SODA3': "/Volumes/DATASETS/SODA3.3.1/OCEAN/",
                    'SODA3_5DAY': "/Volumes/DATASETS/SODA2002/",  # "/cluster/projects/nn9297k/SODA3.3.2/",
                    'NORESM': "/cluster/projects/nn9412k/A20/FORCING/RCP85_ocean/",
                    'GLORYS': "../oceanography/copernicus-marine-data/Global/"}[self.ocean_indata_type]
        except KeyError:
            return KeyError

    def define_atmospheric_forcing_path(self):
        return {'ERA5': "/Volumes/DATASETS/ERA5/"}[self.atmos_indata_type]

    def define_atmos_inputdata_varnames(self):

        return {'ERA5': ['swrad', 'lwrad', 'precip', 'u10', 'v10']}[self.atmos_indata_type]

    def __init__(self):
        logging.info('\n--------------------------\n')
        print('Started ' + time.ctime(time.time()))
        os.environ['WRAP_STDERR'] = 'true'

        # EDIT ===================================================================
        # Set show_progress to "False" if you do not want to see the progress
        # indicator for horizontal interpolation.
        self.show_progress = True
        # Set compileAll to True if you want automatic re-compilation of all the
        # fortran files necessary to run model2roms. Options are "gfortran" or "ifort". Edit
        # compile.py to add other Fortran compilers.
        self.compile_all = False
        # Extract time-series of data for given longitude/latitude
        self.extract_stations = False
        # Define a set of longitude/latitude positions with names to extract into
        # station files (using extractStations)
        if self.extract_stations:
            #  stationNames = ['NorthSea', 'Iceland', 'EastandWestGreenland', 'Lofoten', 'Georges Bank']
            #  lonlist = [2.4301, -22.6001, -47.0801, 13.3801, -67.2001]
            #  latlist = [54.5601, 63.7010, 60.4201, 67.5001, 41.6423]

            self.station_names = ["Ytre Utsira", "Indre Utsira", "Lista"]
            self.latlist = [59.316667, 59.316667, 58.016667]
            self.lonlist = [4.800000, 4.983333, 6.533333]
            self.numberofpoints = 4  # Number of points around lat/lon to extract and average as output

        # Create the bry, init, and clim files for a given grid and input data
        self.create_ocean_forcing = True
        # Create atmospheric forcing for the given grid
        self.create_atmos_forcing = False  # currently in beta stages
        # Create a smaller resolution grid based on your original. Decimates every second for
        # each time run
        self.decimate_gridfile = False
        # Write ice values to file (for Arctic regions)
        self.write_ice = True
        # Write biogeochemistry values to file
        self.write_bcg = False
        # ROMS sometimes requires input of ice and ssh, but if you dont have these write files containing zeros to file
        self.set_2d_vars_to_zero = False
        # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
        self.use_esmf = True
        # Apply filter to smooth the 2D fields after interpolation (time consuming but enhances results)
        self.use_filter = True
        # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
        # Using NETCDF4 automatically turns on compression of files (ZLIB)
        self.output_format = 'NETCDF4'
        self.use_zlib = True
        # Frequency of the input data: usually monthly
        self.time_frequency_inputdata = "month"  # Possible options: "month", "hour", "5days"

        # Path to where results files should be stored
        self.outdir = "../oceanography/model2roms_results/"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)

        # IN GRIDTYPES ------------------------------------------------------------------------------
        # Define what grid type you want to interpolate from (input MODEL data)
        # Currently supported options:
        # 1. NORESM, 2. GLORYS, 3. SODA3, 4. SODA3_5DAY
        self.ocean_indata_type = 'GLORYS'
        self.atmos_indata_type = 'ERA5'

        # Define contact info for final NetCDF files
        self.author_name = "Trond Kristiansen"
        self.author_email = "trond.kristiansen (at) niva.no"

        # Define what grid type you wnat to interpolate from: Can be Z for SIGMA for ROMS
        # vertical coordinate system or ZLEVEL. also define the name of the dimensions in the input files.
        # Options:
        # 1. SIGMA (not properly implemented yet), 2. ZLEVEL
        self.ingrid_type = "SIGMA"  # "ZLEVEL"

        # Define the names of the geographical variables in the input files. These may
        # differ depending how the variable is located in a grid (e.g. Arakawa C grid - ROMS). In
        # SODA 3.3.1 the u and v location is defined by xu_ocean,  yu_ocean while temperature is
        # located in xt_ocean, yt_ocean.
        self.grd_type = 'regular'
        self.lon_name = "longitude"
        self.lat_name = "latitude"
        self.depth_name = "depth"
        self.lon_name_u = "longitude"
        self.lat_name_u = "latitude"
        self.lon_name_v = "longitude"
        self.lat_name_v = "latitude"

        if self.ocean_indata_type == 'SODA3_5DAY':
            self.lon_name = "xt_ocean"
            self.lat_name = "yt_ocean"
            self.depth_name = "st_ocean"
            self.lon_name_u = "xu_ocean"
            self.lat_name_u = "yu_ocean"
            self.lon_name_v = "xu_ocean"
            self.lat_name_v = "yu_ocean"
            self.time_object = []

        self.time_name = "time"
        self.realm = "ocean"
        self.fillvaluein = -1.e20

        # OUT GRIDTYPES ------------------------------------------------------------------------------
        # Define what grid type you want to interpolate to
        # Options: This is just the name of your grid used to identify your selection later
        self.outgrid_name = 'Antarctic'  # "ROHO800", "A20"
        self.outgrid_type = "ROMS"

        # Subset input data. If you have global data you may want to seubset these to speed up reading. Make
        # sure that your input data are cartesian (0-360 or -180:180, -90:90)
        self.subset_indata = False
        if self.subset_indata:
            self.subset = self.define_subset_for_indata()

        # Define nmber of output depth levels
        self.nlevels = 40
        # Define the grid stretching properties (leave default if uncertain what to pick)
        self.vstretching = 4
        self.vtransform = 2
        self.theta_s = 7.0
        self.theta_b = 0.1
        self.tcline = 250.0
        self.hc = 250

        # PATH TO FORCING DATA --------------------------------------------------------------------
        # Define the path to the input data
        self.ocean_forcing_path = self.define_ocean_forcing_data_path()
        self.atmospheric_forcing_path = self.define_atmospheric_forcing_path()

        # PATH TO GRID -----------------------------------------------------------------------------
        # Define the path to the grid file
        self.roms_grid_path = self.define_roms_grid_path()

        # Climatology is only monthly and model2roms needs to know this
        self.isclimatology = False

        # DATE AND TIME DETAILS ---------------------------------------------------------
        # Define the period to create forcing for
        self.start_year = 1993
        self.end_year = 1993
        self.start_month = 1
        self.end_month = 5
        self.start_day = 15
        self.end_day = 31

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

        self.global_varnames = self.define_global_varnames()
        self.input_varnames = self.define_input_data_varnames()
        assert (len(self.global_varnames) == len(self.input_varnames)), "Number and order of global variable " \
                                                                        "names must equal input variable names"

        self.abbreviation = self.define_abbreviation()

        self.clim_name, self.init_name, self.bry_name = self.define_output_filenames()

        if self.isclimatology is True:
            self.clim_name = self.abbreviation + '_' + str(self.ocean_indata_type) + '_climatology.nc'
        self.showinfo()

    def create_grd_objects(self):
        # NO EDIT BELOW ==============================================================================================
        if self.compile_all is True:
            import compile
            compile.compileallgfortran()

        if self.create_atmos_forcing or self.create_ocean_forcing:
            if self.use_esmf:
                try:
                    import ESMF
                except ImportError:
                    raise ImportError("Unable to import ESMF")
                logging.info('[M2R_configRunM2R] Starting logfile for ESMF')
                ESMF.Manager(debug=True)

            # Create the grid object for the output grid
            self.grdROMS = grd.Grd("ROMS", self)
            self.grdROMS.nlevels = self.nlevels
            self.grdROMS.vstretching = self.vstretching
            self.grdROMS.vtransform = self.vtransform
            self.grdROMS.theta_s = self.theta_s
            self.grdROMS.theta_b = self.theta_b
            self.grdROMS.tcline = self.tcline
            self.grdROMS.hc = self.hc
            self.grdROMS.lonname = 'lon_rho'
            self.grdROMS.latname = 'lat_rho'

            self.grdROMS.create_object(self, self.roms_grid_path)
            self.grdROMS.getdims()

            # Create the grid object for the input grid
            self.grdMODEL = grd.Grd("FORCINGDATA", self)
            self.grdMODEL.grdType = self.grd_type
            self.grdMODEL.lonName = self.lon_name
            self.grdMODEL.latName = self.lat_name
            self.grdMODEL.depthName = self.depth_name
            self.grdMODEL.fillval = self.fillvaluein
