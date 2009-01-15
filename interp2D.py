
import numpy as np
import os, sys, string, datetime
from griddata import griddata
from netCDF4 import Dataset
import plotData

__author__   = 'Trond Kristiansen and Bjorn Aadlandsvik'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2008, 12, 18)
__version__  = "1.1"
__status__   = "Development"

                        
                      
def doHorInterpolation(var,grdROMS,grdSODA,data,map,time):
    
    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho
  
    
    tx, ty = map.ll2grid(np.asarray(grdSODA.lon), np.asarray(grdSODA.lat))
    
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))
    
    for t in xrange(time):
        print ' --> Interpolating horziontally for time %s with dimensions (%s,%s)'%(t,Lp,Mp)
        for k in xrange(grdSODA.Nlevels):
        
            xlist = []
            ylist = []
            zlist = []
    
            for i in xrange(len(grdSODA.lon[1,:])):
                for j in xrange(len(grdSODA.lat[:,1])):
                    x = tx[j,i]
                    y = ty[j,i]
                    
                    """
                    Find positions in or near the grid
                    """
                    
                    if ( -20 < x < Lp+10 ) and (-20 < y < Mp+10 ):
                        
                        if data[t,k,j,i]>-9.98e+33:
                            
                            xlist.append(x)
                            ylist.append(y)
                            zlist.append(data[t,k,j,i])
                            #print i,j, data[t,k,j,i]
           #else:
                       #     print 'skipping value %s at %s,%s'%(data[t,k,j,i],y,x)
            #print max(xlist),min(xlist),max(ylist),min(ylist)
            #print 'number of data points in grid %s for time %s for depth %s'%(len(zlist),t,int(grdSODA.Nlevels)-1-k)
            
            Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg, masked=False,fill_value=0)

            if var=='temperature':    
                grdROMS.t[t,k,:,:]=Zg
            elif var=='salinity':
                grdROMS.s[t,k,:,:]=Zg
            elif var=='uvel':
                grdROMS.u[t,k,:,:]=Zg
            elif var=='vvel':
                grdROMS.v[t,k,:,:]=Zg
                
           # plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
   
                     
def doHorInterpolationSSH(var,grdROMS,grdSODA,data,map,time):
    
    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho
  
    
    tx, ty = map.ll2grid(np.asarray(grdSODA.lon), np.asarray(grdSODA.lat))
    
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))
    
    for t in xrange(time):
        print ' --> Interpolating horziontally for time %s with dimensions (%s,%s)'%(t,Lp,Mp)
    
        xlist = []
        ylist = []
        zlist = []

        for i in xrange(len(grdSODA.lon[1,:])):
            for j in xrange(len(grdSODA.lat[:,1])):
                x = tx[j,i]
                y = ty[j,i]
            
                if ( -20 < x < Lp+10 ) and (-20 < y < Mp+10 ):
                    
                    if data[t,j,i]>-9.98e+33:
                        
                        xlist.append(x)
                        ylist.append(y)
                        zlist.append(data[t,j,i])
   
        Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg, masked=False,fill_value=0)

        grdROMS.ssh[t,:,:]=Zg
        
        plotData.contourMap(grdROMS,grdSODA,Zg,"1",var)
                