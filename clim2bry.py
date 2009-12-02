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
__modified__ = datetime(2009, 11,18)
__version__  = "1.2"
__status__   = "Development"

def help():
    """
    This script generates boundary (BRY) files from the climatology (CLIM) files. The
    climatology files are created using createForcing option in main.py of soda2roms package.
    
    Since the variables have different lengths in eta and xi directions, the
    clips of data along the East, WEst, and North and South transects will differ in size. The correct
    sizes are defined below, where No is the number of vertical levels (length of s_rho):
    
    Define the sizes of the clips along xi :
    North and South = salt(No,Lpo),
                      temp(No,Lpo),
                      v(No,Lpo),
                      vbar(1,Lpo),
                      zeta(1,Lpo),
                      u(No,Lp)
                      ubar(1,Lp)
                      
    Define the sizes of the clips along eta :
    East and West   = salt(No,Mpo),
                      temp(No,Mpo),
                      v(No,Mp),
                      vbar(1,Mp),
                      zeta(1,Mpo),
                      u(No,Mpo)
                      ubar(1,Mpo)"""
                      
def writeBry(grdROMS,year,bryName,climName):
    
    """See help function for definition of these variables: """
    
    Lpo = grdROMS.Lp
    Mpo = grdROMS.Mp
    Lp=grdROMS.Lp-1
    Mp=grdROMS.Mp-1
    
    """Open the CLIM file"""
    clim    = Dataset(climName,'r')
    """Generate the BRY netcdf4 file that we will use to fill in data"""
    IOBry.createBryFile(grdROMS,bryName)
    """Now open the file we created"""
    f = Dataset(bryName, mode='a', format='NETCDF4', zlib=True)
     
    """Get the time from the clim file"""
    climtime = np.array(clim.variables["ocean_time"][:])
    ntimes = len(climtime)

    """For each time in CLIM file, save clips of boundary data to BRY file"""
    for itime in range(ntimes):
        
        temp        = np.array(clim.variables["temp"][itime,:,:,:])
        salt        = np.array(clim.variables["salt"][itime,:,:,:])
        ssh         = np.array(clim.variables["zeta"][itime,:,:])
        u           = np.array(clim.variables["u"][itime,:,:,:])
        v           = np.array(clim.variables["v"][itime,:,:,:])
        ubar        = np.array(clim.variables["ubar"][itime,:,:])
        vbar        = np.array(clim.variables["vbar"][itime,:,:])
       
        """Cut the boundary sections along Lpo, and Mpo to create the
        North, South, West, and East boundary clips. Write each clip
        to file for each time step"""
        
        """NOTE: For East and West only V current (v,vbar) have size Mp, others have size Mpo"""
        temp_west = np.squeeze(temp[:,:,0])
        salt_west = np.squeeze(salt[:,:,0])
        ssh_west  = np.squeeze(ssh[:,0])
        u_west    = np.squeeze(u[:,:,0])
        v_west    = np.squeeze(v[:,:,0])
        ubar_west = np.squeeze(ubar[:,0])
        vbar_west = np.squeeze(vbar[:,0])
      
        temp_east = np.squeeze(temp[:,:,Lp])
        salt_east = np.squeeze(salt[:,:,Lp])
        ssh_east  = np.squeeze(ssh[:,Lp])
        u_east    = np.squeeze(u[:,:,Lp-1])
        v_east    = np.squeeze(v[:,:,Lp])
        ubar_east = np.squeeze(ubar[:,Lp-1])
        vbar_east = np.squeeze(vbar[:,Lp])
        
        """NOTE: For South and North only U current (u,ubar) have size Lp, others have size Lpo"""
        temp_south = np.squeeze(temp[:,0,:])
        salt_south = np.squeeze(salt[:,0,:])
        ssh_south  = np.squeeze(ssh[0,:])
        u_south    = np.squeeze(u[:,0,:])
        v_south    = np.squeeze(v[:,0,:])
        ubar_south = np.squeeze(ubar[0,:])
        vbar_south = np.squeeze(vbar[0,:])
        
        temp_north = np.squeeze(temp[:,Mp,:])
        salt_north = np.squeeze(salt[:,Mp,:])
        ssh_north  = np.squeeze(ssh[Mp,:])
        u_north    = np.squeeze(u[:,Mp,:])
        v_north    = np.squeeze(v[:,Mp-1,:])
        ubar_north = np.squeeze(ubar[Mp,:])
        vbar_north = np.squeeze(vbar[Mp-1,:])
         
       
        """------- Write time to file -------------------------------"""
        d= num2date(climtime[itime],units=clim.variables['ocean_time'].long_name,calendar=clim.variables['ocean_time'].calendar)
        print 'clim2bry.py => Appending data to file %s for time %s'%(bryName,d)
            
        f.variables['ocean_time'][itime]      = climtime[itime] 
        
        """------- Write out western boundary variables ------------"""
       
        f.variables['temp_west'][itime,:,:] = temp_west
        f.variables['salt_west'][itime,:,:] = salt_west
        f.variables['zeta_west'][itime,:]   = ssh_west
        
        f.variables['u_west'][itime,:,:]    = u_west
        f.variables['v_west'][itime,:,:]    = v_west
        f.variables['ubar_west'][itime,:]   = ubar_west
        f.variables['vbar_west'][itime,:]   = vbar_west
    
        """------- Write out eastern boundary variables ------------"""
        
        f.variables['temp_east'][itime,:,:] = temp_east
        f.variables['salt_east'][itime,:,:] = salt_east
        f.variables['zeta_east'][itime,:]   = ssh_east
        f.variables['u_east'][itime,:,:]    = u_east
        f.variables['v_east'][itime,:,:]    = v_east
        f.variables['ubar_east'][itime,:]   = ubar_east
        f.variables['vbar_east'][itime,:]   = vbar_east
        
        """------- Write out southern boundary variables ------------"""
    
        f.variables['temp_south'][itime,:,:] = temp_south
        f.variables['salt_south'][itime,:,:] = salt_south
        f.variables['zeta_south'][itime,:]   = ssh_south
        f.variables['u_south'][itime,:,:]    = u_south
        f.variables['v_south'][itime,:,:]    = v_south
        f.variables['ubar_south'][itime,:]   = ubar_south
        f.variables['vbar_south'][itime,:]   = vbar_south
        
        """------- Write out northern boundary variables ------------"""
    
        f.variables['temp_north'][itime,:,:] = temp_north
        f.variables['salt_north'][itime,:,:] = salt_north
        f.variables['zeta_north'][itime,:]   = ssh_north
        f.variables['u_north'][itime,:,:]    = u_north
        f.variables['v_north'][itime,:,:]    = v_north
        f.variables['ubar_north'][itime,:]   = ubar_north
        f.variables['vbar_north'][itime,:]   = vbar_north
        
    f.close()
