<h1>Model2roms</h1>

<a href="https://badge.fury.io/gh/trondkr%2Fmodel2roms"><img src="https://badge.fury.io/gh/trondkr%2Fmodel2roms.svg" alt="GitHub version" height="18"></a>

<a href="https://codeclimate.com/github/codeclimate/codeclimate/maintainability"><img src="https://api.codeclimate.com/v1/badges/a99a88d28ad37a79dbf6/maintainability" /></a>

Model2roms is a Python toolbox for creating the necessary climatology, boundary, and initial forcing files 
required to run the ROMS (<a href="http://myroms.org/" target="_blank">Regional Ocean Modeling System</a>) model. The latest version of model2roms can convert several popular model hindcasts and projections including the NORESM (Norways Earth System Model), SODA global re-analysis, HYCOM, World Ocean Atlas (WOA), and GLORYS (Mercator Ocean) to a use as forcing files for a given ROMS grid structure.

<h3>Latest updates</h3>
<ul>
<li>29.04.2018: Updated support for ESMPy v7.1.0r installed with Anaconda</li>
<li>Support for using Earth System Modeling Framework as the default interpolation method. This allows the input data to be on any kind of grid structure (e.g. irregular) as long as geographical information such as longitude and latitude of grid cells are available. The implementation uses the Python interface to ESMF which can be found here: <a href="https://www.earthsystemcog.org/projects/esmpy/" target="_blank">www.earthsystemcog.org/projects/esmpy/</a>. Using ESMF significantly increases the speed of the interpolation. As an example, interpolating one variable (e.g. temperature distribution) from a global irregular grid to a local non-rectangular grid of size 1250x789, at 70 different depth levels, takes 3 seconds on a Mac Laptop Pro. Additional information as to how to install ESMF and ESMPy on Mac OSX is available <a href="http://www.trondkristiansen.com/?page_id=1302" target="_blank">www.trondkristiansen.com/</a></li>
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
python main.py
```
<h3>Options for interpolation</h3>
The latest version of model2roms has adapated the use of the ESMF python package to handle all of the horizontal interpolations. This has significantly speeded up the interpolations and also solved a number of problems if the grid covers the Nort or South Poles. In addition, ESMF can handle any input type grid and therefore making it very easy to convert any type of model into forcing files for ROMS. However, often the target grid has higher resolution than the source grid which opens up areas (e.g. along the coastlines) where you have no data. Model2roms contains an option `useFilter` that will extrapolate data to fill these areas with no data using a Laplace operator. This is quite useful, but also time-consuming and should be turned off unless you need it:

Without filter            | With filter
:-------------------------:|:-------------------------:
<img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withoutfilter_time_75190.0.png" width=100%>  |  <img src="http://www.trondkristiansen.com/wp-content/gallery/romstools/temperature_depth_ESMF_0_withfilter_time_75190.0.png" width=100%>

<h3>Optional settings</h3>
Prior to run model2roms you have to specify a number of settings so that the program can identify where input and grid files can be found. In addition, you can specify what sort of run you are doing by turning options on and off. Most of the general settings are found in `main.py`, a few definitions for variable names are found in `model2roms.py`, and finally some settings for the grid specifications re found in `grd.py`. Eventually, all of the settings will be moved to one file. Still, the main settings are the following:
``` python
    # Set showprogress to "False" if you do not want to see the progress
    # indicator for horizontal interpolation.
    showprogress = True
    
    # Set compileAll to True if you want automatic re-compilation of all the
    # fortran files necessary to run soda2roms. You need to edit compile.py for this
    compileAll = False

    # Extract time-series of data for given set of longitude/latitudes. Useful to create time-series of input files.
    extractStations = False
    # Define a set of longitude/latitude positions with names to extract into
    
    # Create the bry, init, and clim files for a given ROMS grid and forcing data (e.g. SODA, HYCOM)
    createForcing = True
    # Create a smaller resolution grid based on your original. Decimates every second for
    # each time run. This option is used for large grids where you want to do a number of testing prior ro creating the final      # forcing files.
    decimateGrid = False
    
    # Write ice values to file (for Arctic regions)
    writeIce = True
    
    # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
    useESMF = True
    
    # Apply filter to smooth the 2D fields after interpolation (time consuming)
    useFilter = False

    # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
    # Using NETCDF4 automatically turns on compression of files (ZLIB)
    myformat='NETCDF4'
    
    # Set the input data MODEL type (SODA, SODAMONTHLY,GLORYS2V1,WOAMONTHLY,NORESM)
    mytype = 'NORESM'

    # Define what grid type you wnat to interpolate to (see grd.py for details):
    gridtype = "NS8KM"

    # Define the paths to the input data
    if mytype == 'NORESM':
        modelpath = "/Users/trondkr/Projects/RegScen/NRCP45AERCN_f19_g16_CLE_01/"

    # Define the path to the grid file
    if gridtype == "NS8KM":
        romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
  
    # Define the start year and month and end year and month.
    start_year  = 2017
    end_year    = 2026
    start_month = 1
    end_month   = 12
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
