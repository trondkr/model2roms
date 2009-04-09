
import numpy as np
import datetime
from griddata import griddata
import plotData
import mpl_toolkits.basemap as mp
import geoProjection
import sys
import cl

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2008, 12, 18)
__modified__ = datetime.datetime(2009, 3, 25)
__version__  = "1.3"
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
        print 'Interpolation level %s'%(k)

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
           
        plotData.contourMap(grdROMS,grdSODA,Zg,k,var)
                  
                      
def doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data):
    
    #map = mp.Basemap(resolution='c',projection='robin',lon_0=0.0)
    #map = mp.Basemap(projection='ortho',lat_0=90,lon_0=0.0,
    #              resolution='c',area_thresh=10000.)

    for k in xrange(grdSODA.Nlevels):
       # print 'Interpolating level',k
        
        i0 = np.argmin(np.fabs(grdSODA.lon[0,:]-180))

        dataout = np.zeros(np.squeeze(data[k,:,:]).shape,data.dtype)
        lonsout = np.zeros((len(grdSODA.lon[0,:])),grdSODA.lon.dtype)
        
        # Extract 359 is specific for SODA data as they only have maximum 359  
        lonsout[0:len(grdSODA.lon[0,:])-i0] = grdSODA.lon[0,i0:]-359
       
        lonsout[len(grdSODA.lon[0,:])-i0:] = grdSODA.lon[0,1:i0+1]
        
        dataout[:,0:len(grdSODA.lon[0,:])-i0]  = data[k,:,i0:]
        dataout[:,len(grdSODA.lon[0,:])-i0:] = data[k,:,1:i0+1]
        #
    
        #lons,lats=np.meshgrid(lonsout,np.squeeze(grdSODA.lat[:,0]))
        
        #xin, yin   = map(lons,lats)
        #xout, yout = map(grdROMS.lon_rho,grdROMS.lat_rho)
        
        #xin=xin[0,:]
        #yin=yin[:,0]
      
        Zg = mp.interp(dataout,lonsout,grdSODA.lat[:,0],grdROMS.lon_rho,grdROMS.lat_rho,
                       checkbounds=False, masked=False, order=1)
   
        Zin = np.zeros((grdROMS.lon_rho.shape),dtype=np.float64, order='Fortran')
        Zin = Zg #
     
     #   print '--->cleanArray : Using %s number of points to fill in gaps of data'%(grdROMS.maxDistHorisontal)
        Zin = cl.cleanarray.sweep(np.asarray(grdROMS.depth,order='Fortran'),
                                  float(grdSODA.z_r[k]),
                                  int(grdROMS.minDistPoints),
                                  int(grdROMS.maxval),
                                  int(grdROMS.maxDistHorisontal),
                                  int(grdROMS.maxDistVertical),
                                  np.asarray(Zg,order='Fortran'),
                                  np.asarray(Zin,order='Fortran'),
                                  np.asarray(grdROMS.mask_rho,order='Fortran'),
                                  int(grdROMS.xi_rho),
                                  int(grdROMS.eta_rho))
        
        if var=='temperature':    
            grdROMS.t[k,:,:]=Zin*grdROMS.mask_rho
        elif var=='salinity':
            grdROMS.s[k,:,:]=Zin*grdROMS.mask_rho
        elif var=='uvel':
            grdROMS.u[k,:,:]=Zin*grdROMS.mask_rho
        elif var=='vvel':
            grdROMS.v[k,:,:]=Zin*grdROMS.mask_rho
        
        #plotData.contourMap(grdROMS,grdSODA,Zin*grdROMS.mask_rho,k,var)
      
   

def doHorInterpolationSSHRegularGrid(var,grdROMS,grdSODA,data):

    i0 = np.argmin(np.fabs(grdSODA.lon[0,:]-180))

    dataout = np.zeros(np.squeeze(data[:,:]).shape,data.dtype)
    lonsout = np.zeros((len(grdSODA.lon[0,:])),grdSODA.lon.dtype)
    
    # Extract 359 is specific for SODA data as they only have maximum 359  
    lonsout[0:len(grdSODA.lon[0,:])-i0] = grdSODA.lon[0,i0:]-359
   
    lonsout[len(grdSODA.lon[0,:])-i0:] = grdSODA.lon[0,1:i0+1]
    
    dataout[:,0:len(grdSODA.lon[0,:])-i0]  = data[:,i0:]
    dataout[:,len(grdSODA.lon[0,:])-i0:] = data[:,1:i0+1]
    
    Zg = mp.interp(dataout,lonsout,grdSODA.lat[:,0],grdROMS.lon_rho,grdROMS.lat_rho,
                   checkbounds=False, masked=False, order=1)

    Zin = np.zeros((grdROMS.lon_rho.shape),dtype=np.float64, order='Fortran')
    Zin = Zg
  
    #print '--->cleanArray : Using %s number of points to fill in gaps of data'%(grdROMS.smoothradius)
    Zin = cl.cleanarray.sweep(np.asarray(grdROMS.depth,order='Fortran'),
                                  float(grdSODA.z_r[0]),
                                  int(grdROMS.minDistPoints),
                                  int(grdROMS.maxval),
                                  int(grdROMS.maxDistHorisontal),
                                  int(grdROMS.maxDistVertical),
                                  np.asarray(Zg,order='Fortran'),
                                  np.asarray(Zin,order='Fortran'),
                                  np.asarray(grdROMS.mask_rho,order='Fortran'),
                                  int(grdROMS.xi_rho),
                                  int(grdROMS.eta_rho))
        
   
    grdROMS.ssh[:,:]=Zin*grdROMS.mask_rho
    #plotData.contourMap(grdROMS,grdSODA,Zin,1,var)
        
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
            