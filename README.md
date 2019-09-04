<h1>Model2roms</h1>
<img alt="GitHub" src="https://img.shields.io/github/license/trondkr/model2roms.svg"><img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/trondkr/model2roms.svg">

Model2roms is a Python toolbox for creating the necessary climatology, boundary, and initial forcing files 
required to run the ROMS (<a href="http://myroms.org/" target="_blank">Regional Ocean Modeling System</a>) model. The latest version of model2roms can convert several popular model hindcasts and projections including the NORESM (Norways Earth System Model), SODA global re-analysis, HYCOM, World Ocean Atlas (WOA), and GLORYS (Mercator Ocean) to a use as forcing files for a given ROMS grid structure.

<h3>Latest updates</h3>
<ul>
<li><b>03.09.2019</b>:Added option to use <a href="http://dsrs.atmos.umd.edu/DATA/soda3.3.2/REGRIDED/ocean/" target="_blank"> SODA3 5 day averages</a> as forcing files (<b>SODA3_5DAY</b>). This provides a great opportunity to generate forcing that contains much higher temporal resolution compared to the usual monthly forcing fields. To use this option, define <ul><li>self.indatatype = 'SODA3_5DAY'</li>
        <li>self.timefrequencyofinputdata = "5days"</li>
        </ul>
</li>
        
<li><b>30.04.2019</b>:Fixed progressbar options for Python 3.7</li>
<li><b>15.04.2019</b>:Added `LDFLAG` option in `compile.py` for compiling Fortran files using f2py on Python3.7</li>
<li><b>13.03.2019</b>: Model2roms has been refactored and the code is now easier to read. Added support for interpolating and writing BGC (biogeochemistry) data to the BRY, CLIM, and INIT files. The BGC currently only supports NorESM BGC data as input and the output format (variable names, units) is intended to be used as input to the ERSEM model.</li>
<li>05.09.2018: Model2roms has been refactored and improved in the following way: we now use an object to store all of the configurations (configM2R.py), the run script has been improved (runM2R.py), supports SODA3 and GLORYS2V4 as forcing inputfiles, and severeal minor bugs has been fixed.</li>
<li>Support for using Earth System Modeling Framework as the default interpolation method. This allows the input data to be on any kind of grid structure (e.g. irregular) as long as geographical information such as longitude and latitude of grid cells are available. The implementation uses the Python interface to ESMF which can be found here: <a href="https://www.earthsystemcog.org/projects/esmpy/" target="_blank">www.earthsystemcog.org/projects/esmpy/</a>. Using ESMF significantly increases the speed of the interpolation. As an example, interpolating one variable (e.g. temperature distribution) from a global irregular grid to a local non-rectangular grid of size 1250x789, at 70 different depth levels, takes 3 seconds on a Mac Laptop Pro. Additional information as to how to install ESMF and ESMPy from source on Mac OSX is available <a href="http://www.trondkristiansen.com/?page_id=1302" target="_blank">www.trondkristiansen.com/</a>. For most people, installing using Anaconda would be the best option:

```bash
conda create -n model2roms
conda config --add channels conda-forge
conda activate model2roms
conda install esmpy netcdf4 progressbar
```

</li>

<li>Added support for ICE variables. The latest version writes to file (init, bry, and clim) all necessary ice variables required to run ROMS with ice.</li>
<li>Added support to generate forcing using global <a href="http://sextant.ifremer.fr/record/7a7b31cb-9e7b-4b5b-9ff3-c1165b51f79b/" target="_blank">GLORYS2V3</a> files including sea ice</li>
<li>Many general improvements to the code such as more generic time and date methods.</li>
</ul>

<h2>Background</h2>

For the last few years I have been working on this toolbox model2roms. This toolbox is a collection of Python
and Fortran modules that can be used to create climatology (CLIM), initial (INIT), and boundary (BRY) files 
necessary to run the <a href="www.myroms.org">ROMS</a> model. Currently, the program takes rectangular gridded
forcing files at Z-levels as input. These data are first interpolated to the ROMS grid at Z-levels.
Next the Z-levels are interpolated vertically to the sigma layers
(<a href="https://www.myroms.org/wiki/index.php/Vertical_S-coordinate">S-coordinates</a>).
For U and V velocities, the interpolation is done at RHO points, and then
interpolated to U and V points ((eta_u,xi_u), and (eta_v, xi_v)).
All interpolated values are written to netCDF4 files using compression (zlib) to minimize file size. The result of
running model2roms is one file for each CLIM, INIT, and BRY files.
UBAR and VBAR (barotropic flow) are calculated from U and V velocities. Time is stored as julian
day from 01/01/1948 (see model2roms.py). Make sure to edit the main.py file before you run the toolbox using:

```html
python runM2R.py
```
<h3>Options for interpolation</h3>
The latest version of model2roms has adapated the use of the ESMF python package to handle all of the horizontal interpolations. This has significantly speeded up the interpolations and also solved a number of problems if the grid covers the Nort or South Poles. In addition, ESMF can handle any input type grid and therefore making it very easy to convert any type of model into forcing files for ROMS. However, often the target grid has higher resolution than the source grid which opens up areas (e.g. along the coastlines) where you have no data. Model2roms contains an option `useFilter` that will extrapolate data to fill these areas with no data using a Laplace operator. This is quite useful, but also time-consuming and should be turned off unless you need it:

Without filter            | With filter
:-------------------------:|:-------------------------:
<img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withoutfilter_time_75190.0.png" width=100%>  |  <img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withfilter_time_75190.0.png" width=100%>

<h3>Optional settings</h3>
Prior to run model2roms you have to specify a number of settings so that the program can identify where input and grid files can be found. In addition, you can specify what sort of run you are doing by turning options on and off. All of the user settings are done in `configM2R.py`, a few definitions for variable names are found in `model2roms.py`, and finally a few settings for the grid specifications are found in `grd.py`. Eventually, all of the settings will be moved to one file. Still, the main settings are the following:

```Python

        # Set showprogress to "False" if you do not want to see the progress

        # indicator for horizontal interpolation.
        self.showprogress = False
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
        self.writeice = False
        # Write biogeochemistry values to file
        self.writebcg = False
        # ROMS sometimes requires input of ice and ssh, but if you dont have these write zero files to file
        self.set2DvarsToZero=False
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
        self.indatatype = 'SODA3'

        # Define contact info for final NetCDF files
        self.authorname = "Trond Kristiansen"
        self.authoremail = "trond.kristiansen (at) niva.no"

        # Define what grid type you wnat to interpolate from: Can be Z for SIGMA for ROMS
        # vertical coordinate system or ZLEVEL. also define the name of the dimensions in the input files.
        # Options:
        # 1. SIGMA (not prpoerly implemented yet), 2. ZLEVEL
        self.ingridtype = "SIGMA"

        # Define the names of the geographical variables in the input files
        self.grdtype = 'regular'
        self.lonname = "longitude"
        self.latname = "latitude"
        self.depthname = "depth"
        self.timename = "time"
        self.realm = "ocean"
        self.fillvaluein = -1.e20

        # OUT GRIDTYPES ------------------------------------------------------------------------------
        # Define what grid type you wnat to interpolate to
        # Options: This is just the name of your grid used to identify your selection later
        self.outgrid = "ROHO800"
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
        self.start_year = 1980
        self.end_year = 2014
        self.start_month = 1
        self.end_month = 12
        self.start_day = 15
        self.end_day = 15
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

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
