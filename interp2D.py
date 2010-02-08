# -*- coding: utf-8 -*-
import numpy as np
import datetime
#from griddata import griddata
import plotData
import mpl_toolkits.basemap as mp
import geoProjection
import sys
import cl
import time

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2008, 12, 18)
__modified__ = datetime.datetime(2009, 3, 25)
__version__  = "1.3"
__status__   = "Development"

     
def doHorInterpolationIrregularGrid(var,grdROMS,grdMODEL,data):

    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho

    map=geoProjection.stereographic_wedge(-65.0,52.0,-71.0,47.2,0.15)
         

    tx, ty = map.ll2grid(np.asarray(grdMODEL.lon), np.asarray(grdMODEL.lat))
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))
    
    for k in xrange(grdMODEL.Nlevels):
        print 'Interpolation level %s'%(k)

        xlist = []
        ylist = []
        zlist = []

        for i in xrange(len(grdMODEL.lon[1,:])):
            for j in xrange(len(grdMODEL.lat[:,1])):
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
           
      #  plotData.contourMap(grdROMS,grdMODEL,Zg,k,var)
                  
                      
def doHorInterpolationRegularGrid(var,grdROMS,grdMODEL,data,show_progress):
    
    if show_progress is True:
        from progressBar import progressBar
        # find unicode characters here: http://en.wikipedia.org/wiki/List_of_Unicode_characters#Block_elements
        empty  =u'\u25FD'
        filled =u'\u25FE'

        progress = progressBar(color='red',width=24, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))

    indexROMS_Z_ST = (grdMODEL.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
    array1=np.zeros((indexROMS_Z_ST),dtype=np.float64)

    for k in xrange(grdMODEL.Nlevels):
        
        i0 = np.argmin(np.fabs(grdMODEL.lon[0,:]-180))

        dataout = np.zeros(np.squeeze(data[k,:,:]).shape,data.dtype)
        lonsout = np.zeros((len(grdMODEL.lon[0,:])),grdMODEL.lon.dtype)
        
        lonsout[0:len(grdMODEL.lon[0,:])-i0] = grdMODEL.lon[0,i0:]-360
       
        lonsout[len(grdMODEL.lon[0,:])-i0:] = grdMODEL.lon[0,1:i0+1]

        dataout[:,0:len(grdMODEL.lon[0,:])-i0]  = data[k,:,i0:]
        dataout[:,len(grdMODEL.lon[0,:])-i0:] = data[k,:,1:i0+1]
       
        Zg = mp.interp(dataout,lonsout,grdMODEL.lat[:,0],grdROMS.lon_rho,grdROMS.lat_rho,
                       checkbounds=False, masked=False, order=1)
       
        Zin = np.zeros((grdROMS.lon_rho.shape),dtype=np.float64, order='Fortran')
        Zin = Zg
        Zin = cl.cleanarray.sweep(np.asarray(grdROMS.depth,order='Fortran'),
                                  float(grdMODEL.z_r[k]),
                                  int(grdROMS.minDistPoints),
                                  int(grdROMS.maxval),
                                  int(grdROMS.maxDistHorisontal),
                                  int(grdROMS.maxDistVertical),
                                  np.asarray(Zg,order='Fortran'),
                                  np.asarray(Zin,order='Fortran'),
                                  np.asarray(grdROMS.mask_rho,order='Fortran'),
                                  int(grdROMS.xi_rho),
                                  int(grdROMS.eta_rho))

        Zin=Zin*grdROMS.mask_rho
        
        array1[k,:,:]=Zin
        
        #plotData.contourMap(grdROMS,grdMODEL,Zin,k,var)
        
        if show_progress is True:
            p=int( ((k+1*1.0)/(1.0*grdMODEL.Nlevels))*100.)
        
            progress.render(p)
    
    return array1
    
    

def doHorInterpolationSSHRegularGrid(var,grdROMS,grdMODEL,data):
    
    indexROMS_Z_ST = (grdMODEL.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
    array1=np.zeros((indexROMS_Z_ST),dtype=np.float64)
    
    i0 = np.argmin(np.fabs(grdMODEL.lon[0,:]-180))

    dataout = np.zeros(np.squeeze(data[:,:]).shape,data.dtype)
    lonsout = np.zeros((len(grdMODEL.lon[0,:])),grdMODEL.lon.dtype)
    
    lonsout[0:len(grdMODEL.lon[0,:])-i0] = grdMODEL.lon[0,i0:]-360
   
    lonsout[len(grdMODEL.lon[0,:])-i0:] = grdMODEL.lon[0,1:i0+1]
    
    dataout[:,0:len(grdMODEL.lon[0,:])-i0]  = data[:,i0:]
    dataout[:,len(grdMODEL.lon[0,:])-i0:] = data[:,1:i0+1]
    
    Zg = mp.interp(dataout,lonsout,grdMODEL.lat[:,0],grdROMS.lon_rho,grdROMS.lat_rho,
                   checkbounds=False, masked=False, order=1)

    Zin = np.zeros((grdROMS.lon_rho.shape),dtype=np.float64, order='Fortran')
    Zin = Zg
  
    Zin = cl.cleanarray.sweep(np.asarray(grdROMS.depth,order='Fortran'),
                                  float(grdMODEL.z_r[0]),
                                  int(grdROMS.minDistPoints),
                                  int(grdROMS.maxval),
                                  int(grdROMS.maxDistHorisontal),
                                  int(grdROMS.maxDistVertical),
                                  np.asarray(Zg,order='Fortran'),
                                  np.asarray(Zin,order='Fortran'),
                                  np.asarray(grdROMS.mask_rho,order='Fortran'),
                                  int(grdROMS.xi_rho),
                                  int(grdROMS.eta_rho))
        
   
    array1[0,:,:]=Zin*grdROMS.mask_rho
    #plotData.contourMap(grdROMS,grdMODEL,np.squeeze(array1[0,:,:]),1,var)
    
    return array1

def doHorInterpolationSSHIrregularGrid(var,grdROMS,grdMODEL,data):
    
    Lp=grdROMS.xi_rho
    Mp=grdROMS.eta_rho
  
    tx, ty = map.ll2grid(np.asarray(grdMODEL.lon), np.asarray(grdMODEL.lat))
    
    tx = np.asarray(tx)
    ty = np.asarray(ty)

    Xg, Yg = np.meshgrid(np.arange(Lp), np.arange(Mp))

    xlist = []
    ylist = []
    zlist = []

    for i in xrange(len(grdMODEL.lon[1,:])):
        for j in xrange(len(grdMODEL.lat[:,1])):
            x = tx[j,i]
            y = ty[j,i]
        
            if ( -20 < x < Lp+10 ) and (-20 < y < Mp+10 ):
                
                if data[j,i]>-9.98e+33:
                    
                    xlist.append(x)
                    ylist.append(y)
                    zlist.append(data[j,i])

    Zg = griddata(np.array(xlist), np.array(ylist), np.array(zlist), Xg, Yg)

    grdROMS.ssh[:,:]=Zg
 
    #plotData.contourMap(grdROMS,grdMODEL,Zg,"1",var)
            
