<h1>Model2roms</h1>

This project is a Python toolbox for creating the necessary climatology, boundary, and initial files 
required to run the ROMS (Regional Ocean Modeling System) model. As a first version this toolkit focus on 
converting the SODA, HYCOM, or GLORYS (Mercator Ocean) model runs to a set of forcing files for a given ROMS 
grid structure.

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