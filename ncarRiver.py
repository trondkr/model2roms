import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import numpy as np
import mpl_toolkits.basemap as mp
import plotData

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009,10,29)
__modified__ = datetime(2009,10,29)
__version__  = "1.0"
__status__   = "Development"

def writeRunoff(grdROMS,outfilename,corefile,timeRiver,runoffRiver):
        
    if os.path.exists(outfilename):
        os.remove(outfilename)
        
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a river runoff forcing file for ROMS"
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
    
    v_time = f1.createVariable('runoff_time', 'd', ('time',),zlib=True)
    v_time.long_name = 'Day of year'
    v_time.units = 'day'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time.cycle_length='365.0'
    v_time[:] = timeRiver
  
    v_river = f1.createVariable('runoff', 'f', ('time','eta_rho','xi_rho'),zlib=True)
    v_river.long_name = 'water flux: runoff'
    v_river.units = 'kg meter-2 second-1'
    v_river.time = 'runoff_time';
    v_river[:,:,:] = runoffRiver
    
    f1.close()



def runoff(grdROMS,corepath):
    
    infile='runoff.2004_08_03.nc'
    precipfile = 'ncar_precip.2004_08_03.nc';

    outfile='riverRunoff_core.nc'
    print 'Outputfile is %s'%(outfile)
    if os.path.exists(outfile): os.remove(outfile)

    coreRiverFile = Dataset(corepath+infile,'r')
    coreTimeFile = Dataset(corepath+precipfile,'r')
    
    runoffCORE = np.array(coreRiverFile.variables["Foxx_o_roff"][0,:,:])
    timeCORE   = np.array(coreTimeFile.variables["TIME"][:])
    
    np.ma.masked_greater(runoffCORE>1.e+10,runoffCORE)
    print runoffCORE
    lonCORE = np.arange(0.5,361.5,1)  
    latCORE = np.arange(-89.5,90.5,1)
    """To be consistent with format of nput longitude, all values must be positive"""
    lon=grdROMS.lon_rho
    lat=grdROMS.lat_rho
    lon=np.where(lon<0,lon+360,lon)
    
    runoff = mp.interp(runoffCORE,lonCORE,latCORE,lon,lat,
                       checkbounds=False, masked=False, order=1)
       
   
    #plotData.contourMap(grdROMS,grdROMS,runoff,"1",'runoff')   
    index = 12,runoff.shape[0],runoff.shape[1]
    runoffRiver = np.zeros((index), dtype=np.float64)
    timeRiver   = np.zeros((12), dtype=np.float64)
    
    for mon in range(12):
        runoffRiver[mon,:,:] = runoff
        timeRiver[mon] = timeCORE[mon]
        
    writeRunoff(grdROMS,outfile,corepath+infile,timeRiver,runoffRiver)
    
    




