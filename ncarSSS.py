import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import numpy as np
import mpl_toolkits.basemap as mp
import plotData
import cl2D
from griddata import griddata
from mpl_toolkits.basemap import Basemap

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009,10,30)
__modified__ = datetime(2009,11,3)
__version__  = "1.0"
__status__   = "Development"
def help():
    """This routine uses the griddata module from Jeff Withaker to do
    an irregular 2D interpolation of the SSS to the ROMS grid. The
    reason for this use instead of interp was because the NCAR SSS data
    land matrix is too crude and leaves areas along the coast not interpolated
    when using interp. The griddata routine fixes that problem, but need to be installed:
    
    NOTES:    
    Installing griddata from Jefrey Whitaker :
    http://code.google.com/p/griddata-python/downloads/list
    
    Download asource code and cd into source directory
    python setup.py install
    
    Test the installation by running : python test.py
    
    Trond Kristiansen, 3.11.2009 on the plane from San Francisco to Washington, DC"""
    
def writeSSS(grdROMS,outfilename,corefile,timeSSS,SSS):
        
    if os.path.exists(outfilename):
        os.remove(outfilename)
        
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a Sea Surface Salinity forcing file for ROMS"
    f1.history = 'Created ' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 classic created using SODA2ROMS with %s input files'%(grdROMS.type) 
       
    # Define dimensions
    f1.createDimension('xi_rho',  grdROMS.xi_rho)
    f1.createDimension('eta_rho', grdROMS.eta_rho)
    f1.createDimension('time', None)
   
    vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Longitude at RHO points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_rho
    
    vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Latitude at RHO points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_rho
    
    v_time = f1.createVariable('SSS_time', 'd', ('time',),zlib=True)
    v_time.long_name = 'Day of year'
    v_time.units = 'day'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time.cycle_length='365.25'
    v_time[:] = timeSSS
  
    v_sss = f1.createVariable('SSS', 'f', ('time','eta_rho','xi_rho'),zlib=True)
    v_sss.long_name = 'water flux: SSS'
    v_sss.units = 'kg meter-2 second-1'
    v_sss.time = 'SSS_time';
    v_sss[:,:,:] = SSS
    
    f1.close()



def seaSurfaceSalinity(grdROMS,corepath):
    
    infile='PHC2_salx.2004_08_03.nc'

    outfile='SSS_core.nc'
    print 'Outputfile is %s'%(outfile)
    
    coreSSSFile  = Dataset(corepath+infile,'r')
    
    sssCORE = np.array(coreSSSFile.variables["SALT"][:,:,:])
    timeCORE   = np.array(coreSSSFile.variables["time"][:])
    
    lonCORE = np.array(coreSSSFile.variables["lon"][:]) 
    latCORE = np.array(coreSSSFile.variables["lat"][:])
    
    """To be consistent with format of nput longitude, all values must be positive"""
    lon=grdROMS.lon_rho
    lat=grdROMS.lat_rho
    lon=np.where(lon<0.0,lon+360,lon)
    
    index =  sssCORE.shape[0],grdROMS.lon_rho.shape[0],grdROMS.lon_rho.shape[1]
    SSS = np.zeros((index), dtype=np.float64)
    timeSSS   = np.zeros((timeCORE.shape), dtype=np.float64)
    
    for mon in range(1):
        SSS2D=sssCORE[mon,:,:]
        
        #SSS2D=np.ma.masked_where(SSS2D==-99.0,SSS2D)
       
        #SSS[mon,:,:] = mp.interp(SSS2D,lonCORE,latCORE,lon,lat,
        #               checkbounds=False, masked=False, order=1)
        #SSS[mon,:,:]=np.ma.masked_where(SSS[mon,:,:]==-99.0,SSS[mon,:,:])
        
        #print SSS[mon,:,:]
        #sys.exit()
        map = Basemap(llcrnrlon=-78.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
                      resolution='i',projection='tmerc',lon_0=-70,lat_0=0,area_thresh=10.)
        
        lonC, latC=np.meshgrid(lonCORE,latCORE)
        xi, yi = map(lon,lat)
        x, y = map(lonC,latC)
        print xi
        xlist = []; ylist = []; zlist = []
        #SSS = np.ma.masked_less(SSS,0.0)
        for i in xrange(len(x)):
            for j in xrange(len(y)):
                  
                #if ( -10+x.min() < x[i] < x.max()+10 ) and (-10+y.min() < y[j] < y.max()+10 ):
                    
                if SSS2D[j,i]>-10:
                    #print lonCORE[i], latCORE[j]
                    xlist.append(x[i,i])
                    ylist.append(y[i,i])
                    zlist.append(SSS2D[j,i])
        print xi
        Zin = griddata(np.array(xlist), np.array(ylist), np.array(zlist), xi, yi)
        
        timeSSS[mon] = float(timeCORE[mon])*365.25/12.
        #Zin = np.zeros((grdROMS.lon_rho.shape),dtype=np.float64, order='Fortran')
        #Zin = SSS[mon,:,:]
        plotData.contourMap(grdROMS,grdROMS,Zin,"1",'runoff')
        
        
    writeSSS(grdROMS,outfile,corepath+infile,timeSSS,SSS)
    
    





