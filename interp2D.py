
import numpy as np
import datetime
from griddata import griddata
import plotData
import mpl_toolkits.basemap as mp
import geoProjection
import sys

__author__   = 'Trond Kristiansen and Bjorn Aadlandsvik'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2008, 12, 18)
__modified__ = datetime.datetime(2009, 3, 10)
__version__  = "1.1"
__status__   = "Development"

     
def doHorInterpolationIrregularGrid(var,grdROMS,grdSODA,data):

    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho

    map=geoProjection.stereographic_wedge(-65.0,52.0,-71.0,47.2,0.15)
         

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
        
        Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg)

        if var=='temperature':    
            grdROMS.t[k,:,:]=Zg
        elif var=='salinity':
            grdROMS.s[k,:,:]=Zg
        elif var=='uvel':
            grdROMS.u[k,:,:]=Zg
        elif var=='vvel':
            grdROMS.v[k,:,:]=Zg
           
        #plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
                  
                      
def doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data):
    
   
    map = mp.Basemap(llcrnrlon=-80,llcrnrlat=31.5,urcrnrlon=-55,urcrnrlat=45,
            resolution='h',projection='tmerc',lon_0=-65,lat_0=40)
    
    xin, yin   = map(grdSODA.lon,grdSODA.lat)
    xout, yout = map(grdROMS.lon_rho,grdROMS.lat_rho)
    
    xin=xin[0,:]
    yin=yin[:,0]
   
    
    for k in xrange(grdSODA.Nlevels):

        datain=np.squeeze(data[k,:,:]) 
        #datain = np.ma.masked_values(datain,grdROMS.fill_value)
        
        Zg = mp.interp(datain, xin, yin, xout, yout, checkbounds=False, masked=False, order=1)
        
        
        if var=='temperature':    
            grdROMS.t[k,:,:]=Zg
        elif var=='salinity':
            grdROMS.s[k,:,:]=Zg
        elif var=='uvel':
            grdROMS.u[k,:,:]=Zg
        elif var=='vvel':
            grdROMS.v[k,:,:]=Zg
       
        #plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
   

def doHorInterpolationSSHRegularGrid(var,grdROMS,grdSODA,data):

    """
    Use basemap instances to project the latitude/longitude pairs from input grid and output
    grid to a general mercator projection. This enables us to use bilinear/nearest neighbor interpolation
    when the input grid is regular.
    """
    map = mp.Basemap(llcrnrlon=-80,llcrnrlat=31.5,urcrnrlon=-55,urcrnrlat=45,
            resolution='h',projection='tmerc',lon_0=-65,lat_0=40)
    
    xin, yin   = map(grdSODA.lon,grdSODA.lat)
    xout, yout = map(grdROMS.lon_rho,grdROMS.lat_rho)
    
    xin=xin[0,:]
    yin=yin[:,0]
   
    
    print 'interpolating field %s at surface'%(var)
    datain=data[:,:] 
    datain = np.ma.masked_values(datain,grdROMS.fill_value)
    
    Zg = mp.interp(datain, xin, yin, xout, yout, checkbounds=False, masked=True, order=1)
  
    grdROMS.ssh[:,:]=Zg
    #plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
        
def doHorInterpolationSSHIrregularGrid(var,grdROMS,grdSODA,data):
    
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

    Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg)

    grdROMS.ssh[:,:]=Zg
 
    #plotData.contourMap(grdROMS,grdSODA,Zg,"1",var)
            