
import numpy as np
import os, sys, string, datetime
from griddata import griddata
import pylab as p
from netCDF4 import Dataset
import plotData

__author__   = 'Trond Kristiansen and Bjorn Aadlandsvik'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2008, 12, 4)
__version__  = "1.0"
__status__   = "Development"

                        
                      
def doHorInterpolation(grdROMS,data,lon_in,lat_in,Lpo,Mpo,map,time,Nlevels):
    
    tx, ty = map.ll2grid(np.asarray(lon_in), np.asarray(lat_in))
    
    tx = np.asarray(tx)
    ty = np.asarray(ty)
    
    mask=np.where(data==-9.99e+33, 0, 1)
    data=data*mask
    
    Xg, Yg = np.meshgrid(np.arange(Lpo), np.arange(Mpo))
    
    for t in xrange(time):
        for k in xrange(4) : #Nlevels):
        
            xlist = []
            ylist = []
            zlist = []
    
            for i in xrange(len(lon_in[1,:])):
                for j in xrange(len(lat_in[:,1])):
                    x = tx[j,i]
                    y = ty[j,i]
                    
                    """
                    Find positions in or near the grid
                    """
                    
                    if ( -1 < x < Lpo+2 ) and (-1 < y < Mpo+2 ):
                        xlist.append(x)
                        ylist.append(y)
                        zlist.append(data[t,k,j,i])
                        
            print 'number of data points in grid %s for time %s for depth %s'%(len(zlist),t,k)
         
            Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg, masked=False)
        
            grdROMS.t[t,k,:,:]=Zg
      
   # plotData.contourMap(grd,grd.depth)
   