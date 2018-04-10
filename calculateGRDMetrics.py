import time
from datetime import datetime
import model2roms
import IOstation
import clim2bry
import decimateGrid
import grd
import numpy as np
import atmosForcing

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2015, 7, 7)
__modified__ = datetime(2015, 7, 7)
__version__ = "1.0"
__status__ = "Development"

doc = """ dist = greatCircle(lon1,lat1,lon2,lat2]


  greatCircle:  Coumputes great circle distance between to points

  This function computes great circle distance between two longitude and latitude points. 
  The Earth is assumed to be a sphere. This approximation is valid for short distances.

  On Input:  Longitude is positive to the east and negative to the
             west.  Latitude is positive to the north and negative
             to the south.

     lon1     longitude point 1 [decimal degrees]
     lat1     latitude  point 1 [decimal degrees]
     lon2     longitude point 2 [decimal degrees]
     lat2     latitude  point 2 [decimal degrees]

  On Output:

     dist     great circle distance between point 1 and point 2
                [meters]

 Adapted from routine written by Pat J. Haley [Harvard University].

"""

def greatCircle(lon1,lat1,lon2,lat2):

	#----------------------------------------------------------------------------
	# Set often used parameters.
	#----------------------------------------------------------------------------

	radius  = 6371.315
	deg2rad = np.pi/180.0
	rad2deg = 180.0/np.pi
	
	#----------------------------------------------------------------------------
	# Convert to radians.
	#----------------------------------------------------------------------------

	slon = np.multiply(lon1,deg2rad)
	slat = np.multiply(lat1,deg2rad)
	elon = np.multiply(lon2,deg2rad)
	elat = np.multiply(lat2,deg2rad)

	#----------------------------------------------------------------------------
	# Compute distance along great circle [kilometers].
	#----------------------------------------------------------------------------

	alpha = np.multiply(np.sin(slat),np.sin(elat)) + np.cos(slat)*np.cos(elat)*np.cos(elon-slon)

	alpha=np.arccos(alpha)
	#km2meter=1000.
	dist=np.multiply(radius,alpha)#*km2meter

	return dist

def calculateGridMetrics(G, GreatCircle, decimate, startindex, endindex):
	doc2 = """function [pm, pn, dndx, dmde] = grid_metrics[G, GreatCircle]


	  GRID_METRICS:  Compute ROMS Grid horizontal metrics

	  [pm, pn, dndx, dmde] = grid_metrics[G, GreatCircle]

	  decimate is used for createing a smaller grid basedon the input data. 
	  The new grid uses every 'decimate' point to calculate the metrics.

	  This function computes horizontal grid spacing metrics from
	  Grid NetCDF file or Grid structure G.

	  On Input:

	    G            A ROMS grid structure as defined in grd.py

	    GreatCircle  Switch indicating how to compute the grid distance:
	                   GreatCircle = true     Great-circle distance
	                   GreatCircle = false    Cartesian distance

	  On Output:

	    pm         Curvilinear coordinate metric in the XI-direction
	                 [1/meters dx = 1/pm]

	    pm         Curvilinear coordinate metric in the ETA-direction
	                 [1/meters dy = 1/pn]

	    dndx       XI-derivative  of inverse metric factor pn [meters],
	                 d[pn]/d[XI]

	    dmde       ETA-derivative of inverse metric factor pm [meters],
	                 d[pm]/d[ETA]

	  If G is a Grid structure and GreatCircle=true, the following values are
	  needed to compute the great circle distances [G.spherical must be 1]:

	    G.spherical     Grid spherical flag [0 | 1]
	    G.lon_rho       RHO-points longitude [decimal degrees]
	    G.lat_rho       RHO-points latitude  [decimal degrees]
	    G.lon_u         U-points   longitude [decimal degrees]
	    G.lat_u         U-points   latitude  [decimal degrees]
	    G.lon_v         V-points   longitude [decimal degrees]
	    G.lat_v         V-points   latitude  [decimal degrees]

	                    longitude:  positive East,   negative West
	                    latitude:   positive North,  negative South

	 Otherwise, if G is a grid structure and GreatCircle=false, the following
	 values are needed to compute Cartesian distances regardless the value of
	 G.spherical:

	    G.spherical     Grid spherical flag [0 | 1]
	    G.x_rho         RHO-points X-location [meters]
	    G.y_rho         RHO-points Y-location [meters]
	    G.x_u           U-points   X-location [meters]
	    G.x_u           U-points   Y-location [meters]
	    G.x_v           V-points   X-location [meters]
	    G.x_v           V-points   Y-location [meters]

	 =========================================================================%
	  Copyright [c] 2002-2014 The ROMS/TOMS Group                            %
	    Licensed under a MIT/X style license                                 %
	    See License_ROMS.txt                           Hernan G. Arango      %
	 =========================================================================%

	"""
	 # Get Grid coordinates.

	
	spherical = G.spherical

	if (GreatCircle is True and spherical=="T"):
	  Xr = G.lon_rho[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yr = G.lat_rho[startindex:endindex:decimate,startindex:endindex:decimate]
	  Xu = G.lon_u[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yu = G.lat_u[startindex:endindex:decimate,startindex:endindex:decimate]
	  Xv = G.lon_v[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yv = G.lat_v[startindex:endindex:decimate,startindex:endindex:decimate]
	else:
	  Xr = G.x_rho[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yr = G.y_rho[startindex:endindex:decimate,startindex:endindex:decimate]
	  Xu = G.x_u[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yu = G.y_u[startindex:endindex:decimate,startindex:endindex:decimate]
	  Xv = G.x_v[startindex:endindex:decimate,startindex:endindex:decimate]
	  Yv = G.y_v[startindex:endindex:decimate,startindex:endindex:decimate]
	
	#----------------------------------------------------------------------------  
	#  Compute grid spacing [meters].
	#----------------------------------------------------------------------------  
	Mp=len(Xr[:,0])
	Lp=len(Xr[0,:])
	L = Lp; Lm = L-1
	M = Mp; Mm = M-1
	dx = np.zeros((Mp,Lp))
	dy = np.zeros((Mp,Lp))
	# Compute grid spacing.

	if (GreatCircle and spherical):
	 
	  dx[1:M,0:Lp] = greatCircle(Xr[0:Mm,0:Lp], Yr[0:Mm,0:Lp],                
	                         Xr[1:M ,0:Lp], Yr[1:M ,0:Lp])
	   	  
	  dx[0  ,0:Lp] = greatCircle(Xr[0   ,0:Lp], Yr[0   ,0:Lp],                
	                         Xr[1   ,0:Lp], Yr[1   ,0:Lp])
	
	  dx[-1 ,0:Lp] = greatCircle(Xr[-2  ,0:Lp], Yr[-2  ,0:Lp],                
	                         Xr[-1   ,0:Lp], Yr[-1   ,0:Lp]) 
		
	  dy[0:Mp,1:Lp] = greatCircle(Xr[0:Mp,0:Lm], Yr[0:Mp,0:Lm],                
	                         Xr[0:Mp,1:Lp ], Yr[0:Mp,1:Lp ])
 		
	  dy[0:Mp,0  ] = greatCircle(Xr[0:Mp,0   ], Yr[0:Mp,0   ],                
	                         Xr[0:Mp,1   ], Yr[0:Mp,1   ])
	  
	  dy[0:Mp,-1] = greatCircle(Xr[0:Mp,-2  ], Yr[0:Mp,-2  ],                
	                         Xr[0:Mp, -1   ], Yr[0:Mp,-1   ])
	
	  dx = dx * 1000.        #great circle function computes
	  dy = dy * 1000.        #distances in meters
	  
	# Compute inverse grid spacing metrics.
	
	pm = 1.0/dx
	pn = 1.0/dy

	#----------------------------------------------------------------------------  
	# Compute inverse metric derivatives.
	#----------------------------------------------------------------------------  

	dndx = np.zeros(np.shape(Xr))
	dmde = np.zeros(np.shape(Xr))
	
	dndx[1:Mp-1,1:Lp-1] = 0.5*(1.0/pn[2:Mp,1:Lp-1] - 1.0/pn[0:Mp-2,1:Lp-1])
	dmde[1:Mp-1,1:Lp-1] = 0.5*(1.0/pm[1:Mp-1,2:Lp] - 1.0/pm[1:Mp-1,1:Lp-1])
	
	return dndx,dmde,pm,pn

