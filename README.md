<h1>model2roms</h1>

[![codebeat badge](https://codebeat.co/badges/cf0f4bfc-a186-4dfc-9e6a-291b7e985492)](https://codebeat.co/projects/github-com-trondkr-model2roms-master)
[![Build status](https://badge.buildkite.com/8b25e9518aca90534a2a0755dd0623dc91533ba7ca5bd42dde.svg)](https://buildkite.com/rask-dev-llc/model2roms)
![License][image-4]
[![DOI](https://zenodo.org/badge/11505338.svg)](https://zenodo.org/badge/latestdoi/11505338)

[image-1]:    https://buildkite.com/rask-dev-llc/model2roms

[image-2]:    https://codebeat.co/projects/github-com-trondkr-model2roms-master

[image-3]:    https://badge.buildkite.com/8b25e9518aca90534a2a0755dd0623dc91533ba7ca5bd42dde.svg

[image-4]:    https://img.shields.io/github/last-commit/trondkr/model2roms.svg

Model2roms is a Python toolbox for creating the necessary climatology, boundary, and initial forcing files required to
run the ROMS (<a href="http://myroms.org/" target="_blank">Regional Ocean Modeling System</a>) model. The latest version
of model2roms can convert several popular model hindcasts and projections including the NORESM (Norways Earth System
Model), SODA global re-analysis, HYCOM, World Ocean Atlas (WOA), and GLORYS (Mercator Ocean) to a use as forcing files
for a given ROMS grid structure.

<h3>Introduction</h3>

Model2roms is a Python toolbox for creating the necessary climatology, boundary, and initial forcing files required to
run the ROMS (<a href="http://myroms.org/" target="_blank">Regional Ocean Modeling System</a>) model. The latest version
of model2roms can convert several popular model hindcasts and projections including the NORESM
(Norways Earth System Model)
, **[SODA](https://climatedataguide.ucar.edu/climate-data/soda-simple-ocean-data-assimilation)**
global re-analysis, and
**[GLORYS](https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=GLOBAL_REANALYSIS_PHY_001_030)**
(Mercator Ocean) to be used as forcing files for a given ROMS grid structure.

Model2roms uses the Earth System Modeling Framework (ESMF) as the default interpolation method. This allows the input
data to be on any kind of grid structure (e.g. irregular) provided that geographical information such as the longitude
and latitude of grid cells is available. The implementation uses the Python interface to ESMF which can be found
here: <a href="https://www.earthsystemcog.org/projects/esmpy/" target="_blank">www.earthsystemcog.org/projects/esmpy/</a>
. Using ESMF significantly increases the speed of the interpolation. As an example, interpolating one variable
(e.g. temperature distribution) from a global irregular grid to a local non-rectangular grid of size 1250x789, at 70
different depth levels, takes 3 seconds on a Mac Laptop Pro. For most people, installing using Anaconda would be the
best option to install the required packages to run <b>model2roms</b>. The minimum installation of required packages:

```bash
conda create -n model2roms`
conda config --add channels conda-forge
conda activate model2roms
conda install esmf xarray netcdf4 progressbar2 
```

![Sea-ice concentration Antarctica](https://github.com/trondkr/model2roms/blob/master/Examples/Figures/temp_Antarctic.png)

Model2roms has been developed over several years, usually improved every time a new model configuration of ROMS has been
required for my work. But the user groups of model2roms has also increased over these years and several people have
pointed out bugs and ways to improve this package. The toolbox consists of collection of Python and Fortran modules that
can be used to create climatology (CLIM), initial (INIT), and boundary (BRY) files necessary to run
the <a href="www.myroms.org">ROMS</a> model.

Currently, model2roms takes rectangular gridded forcing files at Z-levels as input. These data are first interpolated to
the ROMS grid at Z-levels. Next the Z-levels are interpolated vertically to the sigma layers
(<a href="https://www.myroms.org/wiki/index.php/Vertical_S-coordinate">S-coordinates</a>). For U and V velocities, the
interpolation is done at RHO points, and then interpolated to U and V points ((eta_u,xi_u), and (eta_v, xi_v)). All
interpolated values are written to netCDF4 files using compression (zlib) to minimize file size. The result of running
model2roms is one file for each CLIM, INIT, and BRY files. UBAR and VBAR (barotropic flow) are calculated from U and V
velocities. Time is stored as julian day from 01/01/1948 (see model2roms.py)

<h4>Compile the Fortran modules</h4>
To get started, compile the Fortran functions into modules callable by Python. First edit the file compile.py and select
your Fortran compiler (currently gfortran and Intel Fortran compiler supported):

```python
python
compile.py
```

Make sure that this successfully creates modules (.so files) that Python can import. Some users have reported that they
have to run each compile command individually to compile (this may depend on your machine and OS). The use of Fortran
modules as part of the calculations significantly speeds up the calculations.

<h4>Running model2roms</h4>
Before you run model2roms you have to edit the configuration file `configM2R.py` to correctly point to the path of your
gridfile, the type of forcing you want and the variables to use.

Once everything is correctly setup you can run model2roms with the command:
```python
python
runM2R.py
```

<h4>Options for interpolation</h4>
Model2roms makes use of the ESMF python package to handle all the horizontal interpolations. This has significantly sped
up the time used on interpolation and also made interpolation more robust across the poles.  
In addition, ESMF can handle any input type grid and therefore making it very easy to convert any type of model into
forcing files for ROMS. However, often the target grid has higher resolution than the source grid which opens up areas
(e.g. along the coastlines) where you have no data. Model2roms contains an option `use_filter` that will extrapolate
data to fill these areas with no data using a Laplace operator. This is quite useful, but also time-consuming and should
be turned off unless you need it:

Without filter            | With filter
:-------------------------:|:-------------------------:
<img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withoutfilter_time_75190.0.png" width=100%>  |  <img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withfilter_time_75190.0.png" width=100%>

<h3>Optional settings</h3>
Prior to run model2roms you have to specify a number of settings so that the program can identify where input and grid
files can be found. In addition, you can specify what sort of run you are doing by turning options on and off. All of
the user settings are done in `configM2R.py`, a few definitions for variable names are found in `model2roms.py`, and
finally a few settings for the grid specifications are found in `grd.py`. Eventually, all of the settings will be moved
to one file. Still, the main settings are the following:

```Python
    def __init__(self):
    logging.info('\n--------------------------\n')
    logging.info('Started ' + time.ctime(time.time()))
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
    self.write_ice = False

    # Write biogeochemistry values to file
    self.write_bcg = False

    # ROMS sometimes requires input of ice and ssh, but if you dont have these write files containing zeros to file
    self.set_2d_vars_to_zero = False

    # Apply filter to smooth the 2D fields after interpolation (time consuming but enhances results)
    self.use_filter = True

    # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
    # Using NETCDF4 automatically turns on compression of files (ZLIB)
    self.output_format = 'NETCDF4'
    self.use_zlib = True

    # Frequency of the input data: usually monthly
    self.time_frequency_inputdata = "month"  # Possible options: "month", "hour", "5days"

    # Path to where results files should be stored
    self.outdir = "../oceanography/NAUTILOS/"
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
    self.ingrid_type = "ZLEVEL"  # "ZLEVEL"

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
    self.fillvaluein = -32767

    # OUT GRIDTYPES ------------------------------------------------------------------------------
    # Define what grid type you want to interpolate to
    # Options: This is just the name of your grid used to identify your selection later
    self.outgrid_name = 'ROHO160'  # "ROHO800", "A20"
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
    self.start_year = 2017
    self.end_year = 2019
    self.start_month = 1
    self.end_month = 11
    self.start_day = 15
    self.end_day = 31
```  

<p style="clear: both;">

<h2>Contact</h2>
<ul>
<li>me @ trondkristiansen.com</li>
<li>http://github.com/trondkr</li>
<li>www.trondkristiansen.com</li>
</ul>
Please send a comment on suggestions, questions, or modifications and improvements you would like to
make to me @ trondkristiansen.com. I would very much like to see this project go much further so that it can be
useful for a variety of purposes.

<h2>License</h2>
The MIT License (MIT)

Copyright (c) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

<h2>Contact</h2>
<ul>
<li>me @ trondkristiansen.com</li>
<li>http://github.com/trondkr</li>
<li>www.trondkristiansen.com</li>
</ul>
Please send a comment on suggestions, questions, or modifications and improvements you would like to
make to me @ trondkristiansen.com. I would very much like to see this project go much further so that it can be
useful for a variety of purposes.

<h2>License</h2>
The MIT License (MIT)
Copyright (c) <year> <copyright holders>
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

