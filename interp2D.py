
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

                        
                      
def doHorInterpolation(var,grdROMS,grdSODA,data,map):
    
    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho

    tx, ty = map.ll2grid(np.asarray(grdSODA.lon), np.asarray(grdSODA.lat))
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))
    
    for k in xrange(grdSODA.Nlevels):
       # print 'Interpolation level %s'%(k)
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
                    
                    if data[k,j,i]>-9.98e+33:
                        
                        xlist.append(x)
                        ylist.append(y)
                        zlist.append(data[k,j,i])
                        #print i,j, data[t,k,j,i]
       #else:
                   #     print 'skipping value %s at %s,%s'%(data[t,k,j,i],y,x)
        
        Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg, masked=False,fill_value=0)

        if var=='temperature':    
            grdROMS.t[k,:,:]=Zg
        elif var=='salinity':
            grdROMS.s[k,:,:]=Zg
        elif var=='uvel':
            grdROMS.u[k,:,:]=Zg
        elif var=='vvel':
            grdROMS.v[k,:,:]=Zg
           
        #plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
   
                     
def doHorInterpolationSSH(var,grdROMS,grdSODA,data,map):
    
    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho
  
    tx, ty = map.ll2grid(np.asarray(grdSODA.lon), np.asarray(grdSODA.lat))
    
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))

    xlist = []
    ylist = []
    zlist = []

    for i in xrange(len(grdSODA.lon[1,:])):
        for j in xrange(len(grdSODA.lat[:,1])):
            x = tx[j,i]
            y = ty[j,i]
        
            if ( -20 < x < Lp+10 ) and (-20 < y < Mp+10 ):
                
                if data[j,i]>-9.98e+33:
                    
                    xlist.append(x)
                    ylist.append(y)
                    zlist.append(data[j,i])

    Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg, masked=False,fill_value=0)

    grdROMS.ssh[:,:]=Zg
 
    #plotData.contourMap(grdROMS,grdSODA,Zg,"1",var)
            