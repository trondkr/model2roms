import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import numpy as np
import IOBry

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 3,2)
__modified__ = datetime(2009, 3,2)
__version__  = "1.0"
__status__   = "Development"


def help():
    
    """
    This script generates boundary files from the climatology files
    created using soda2roms.
    """
    
def writeBry(grdROMS,year,climname):
    
    
    # Loop over all clim files in list 
    #filelist = [clim file list ] #strcat('/Users/trond/ROMS/LA_04/GOM/brylst_',year,'.txt');

    No = len(grdROMS.s_rho)
    Lpo = grdROMS.Lp
    Mpo = grdROMS.Mp
    Lp=grdROMS.Lp-1
    Mp=grdROMS.Mp-1
    
    #for file in filelist:
    clim    = Dataset(climname,'r')
    bryname = 'bry_'+str(climname)
    
    """Generate the BRY netcdf4 file that we will use to fill in data"""
    IOBry.createBryFile(grdROMS,bryname)
    """Now open the file we created"""
    f = Dataset(bryname, mode='a', format='NETCDF4', zlib=True)
     
    climtime = np.array(clim.variables["clim_time"][:])
    
    ntimes = len(climtime)

    for itime in range(ntimes):
        
        temp        = np.array(clim.variables["temp"][itime,:,:,:])
        salt        = np.array(clim.variables["salt"][itime,:,:,:])
        ssh         = np.array(clim.variables["zeta"][itime,:,:])
        u           = np.array(clim.variables["u"][itime,:,:,:])
        v           = np.array(clim.variables["v"][itime,:,:,:])
        ubar        = np.array(clim.variables["ubar"][itime,:,:])
        vbar        = np.array(clim.variables["vbar"][itime,:,:])
        
      
        temp_west = temp[:,:,0].reshape(No,Mpo,order='F').copy()
        salt_west = salt[:,:,0].reshape(No,Mpo,order='F').copy()
        ssh_west  = ssh[:,0].reshape(1,Mpo,order='F').copy()
        u_west    = u[:,:,0].reshape(No,Mpo,order='F').copy()
        v_west    = v[:,:,0].reshape(No,Mp,order='F').copy()
        ubar_west = ubar[:,0].reshape(1,Mpo,order='F').copy()
        vbar_west = vbar[:,0].reshape(1,Mp,order='F').copy()
       
        temp_east = temp[:,:,Lp].reshape(No,Mpo,order='F').copy()
        salt_east = salt[:,:,Lp].reshape(No,Mpo,order='F').copy()
        ssh_east  = ssh[:,Lp].reshape(1,Mpo,order='F').copy()
        u_east    = u[:,:,Lp-1].reshape(No,Mpo,order='F').copy()
        v_east    = v[:,:,Lp].reshape(No,Mp,order='F').copy()
        ubar_east = ubar[:,Lp-1].reshape(1,Mpo,order='F').copy()
        vbar_east = vbar[:,Lp].reshape(1,Mp,order='F').copy()
        
        temp_south = temp[:,0,:].reshape(No,Lpo,order='F').copy()
        salt_south = salt[:,0,:].reshape(No,Lpo,order='F').copy()
        ssh_south  = ssh[0,:].reshape(1,Lpo,order='F').copy()
        u_south    = u[:,0,:].reshape(No,Lp,order='F').copy()
        v_south    = v[:,0,:].reshape(No,Lpo,order='F').copy()
        ubar_south = ubar[0,:].reshape(1,Lp,order='F').copy()
        vbar_south = vbar[0,:].reshape(1,Lpo,order='F').copy()
       
        temp_north = temp[:,Mp,:].reshape(No,Lpo,order='F').copy()
        salt_north = salt[:,Mp,:].reshape(No,Lpo,order='F').copy()
        ssh_north  = ssh[Mp,:].reshape(1,Lpo,order='F').copy()
        u_north    = u[:,Mp,:].reshape(No,Lp,order='F').copy()
        v_north    = v[:,Mp-1,:].reshape(No,Lpo,order='F').copy()
        ubar_north = ubar[Mp,:].reshape(1,Lp,order='F').copy()
        vbar_north = vbar[Mp-1,:].reshape(1,Lpo,order='F').copy()
    
        # ------- Write time to file -------------------------------
        d= num2date(climtime[itime],units=clim.variables['clim_time'].long_name,calendar=clim.variables['clim_time'].calendar)
        print 'Appending data to file %s for time %s'%(bryname,d)
            
        f.variables['bry_time'][itime]      = climtime[itime]
        
        
        # ------- Write out western boundary variables ------------
       
        f.variables['temp_west'][itime,:,:]  = temp_west
        f.variables['salt_west'][itime,:,:] = salt_west
        f.variables['zeta_west'][itime,:]   = ssh_west
        f.variables['u_west'][itime,:,:]    = u_west
        f.variables['v_west'][itime,:,:]    = v_west
        f.variables['ubar_west'][itime,:] = ubar_west
        f.variables['vbar_west'][itime,:] = vbar_west
    
    
        # ------- Write out eastern boundary variables ------------
        
        f.variables['temp_east'][itime,:,:]  = temp_east
        f.variables['salt_east'][itime,:,:] = salt_east
        f.variables['zeta_east'][itime,:]   = ssh_east
        f.variables['u_east'][itime,:,:]    = u_east
        f.variables['v_east'][itime,:,:]    = v_east
        f.variables['ubar_east'][itime,:] = ubar_east
        f.variables['vbar_east'][itime,:] = vbar_east
        
        # ------- Write out southern boundary variables ------------
    
        f.variables['temp_south'][itime,:,:]  = temp_south
        f.variables['salt_south'][itime,:,:] = salt_south
        f.variables['zeta_south'][itime,:]   = ssh_south
        f.variables['u_south'][itime,:,:]    = u_south
        f.variables['v_south'][itime,:,:]    = v_south
        f.variables['ubar_south'][itime,:] = ubar_south
        f.variables['vbar_south'][itime,:] = vbar_south
        
         # ------- Write out northern boundary variables ------------
    
        f.variables['temp_north'][itime,:,:]  = temp_north
        f.variables['salt_north'][itime,:,:] = salt_north
        f.variables['zeta_north'][itime,:]   = ssh_north
        f.variables['u_north'][itime,:,:]    = u_north
        f.variables['v_north'][itime,:,:]    = v_north
        f.variables['ubar_north'][itime,:] = ubar_north
        f.variables['vbar_north'][itime,:] = vbar_north
        
    f.close()