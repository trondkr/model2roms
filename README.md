<h1>Model2roms</h1>

Model2roms is a Python toolbox for creating the necessary climatology, boundary, and initial forcing files 
required to run the ROMS (<a href="http://myroms.org/" target="_blank">Regional Ocean Modeling System</a>) model. The latest version of model2roms can convert several popular model hindcasts and projections including the NORESM (Norways Earth System Model), SODA global re-analysis, HYCOM, World Ocean Atlas (WOA), and GLORYS (Mercator Ocean) to a use as forcing files for a given ROMS grid structure.

<h3>Latest updates</h3>
<ul>
<li>Support for using Earth System Modeling Framework as the default interpolation method. This allows the input data to be on any kind of grid structure (e.g. irregular) as long as geographical information such as longitude and latitude of grid cells are available. The implementation uses the Python interface to ESMF which can be found here: <a href="https://www.earthsystemcog.org/projects/esmpy/" target="_blank">www.earthsystemcog.org/projects/esmpy/</a>. Using ESMF significantly increases the speed of the interpolation. As an example, interpolating one variable (e.g. temperature distribution) from a global irregular grid to a local non-rectangular grid of size 1250x789, at 70 different depth levels, takes 3 seconds on a Mac Laptop Pro. Additional information as to how to install ESMF and ESMPy on Mac OSX is available <a href="http://www.trondkristiansen.com/?page_id=1302" target="_blank">www.trondkristiansen.com/</a></li>
<li>Added support for ICE variables. The latest version writes to file (init, bry, and clim) all necessary ice variables required to run ROMS with ice.</li>
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

Model2roms is available under the MIT license. See the LICENSE file for more info.
