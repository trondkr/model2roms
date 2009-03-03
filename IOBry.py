import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 3,2)
__modified__ = datetime(2009, 3,2)
__version__  = "0.1"
__status__   = "Development"


def help ():
     """
     This function generates a BRY file from scratch. The variables are all created
     for East, West, North, and South. Varibales include:
     salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is bry_time which is days
     since 1948/1/1.
     """

def createBryFile(grdROMS,outfilename):
     
     if os.path.exists(outfilename):
         os.remove(outfilename)
     print 'Creating initial Boundary (BRY) file %s'%(outfilename)
     
     f1 = Dataset(outfilename, mode='w', format='NETCDF4')
     f1.description="This is a boundary file for ROMS"
     f1.history = 'Created ' + time.ctime(time.time())
     f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
     f1.type='NetCDF4 classic created using SODA2ROMS' 
        
     # Define dimensions
     f1.createDimension('xi_rho',  grdROMS.xi_rho)
     f1.createDimension('eta_rho', grdROMS.eta_rho)
     f1.createDimension('xi_u',    grdROMS.xi_u)
     f1.createDimension('eta_u',   grdROMS.eta_u)
     f1.createDimension('xi_v',    grdROMS.xi_v)
     f1.createDimension('eta_v',   grdROMS.eta_v)
     f1.createDimension('xi_psi',    grdROMS.xi_psi)
     f1.createDimension('eta_psi',   grdROMS.eta_psi)
     f1.createDimension('bry_time', None)
     f1.createDimension('s_rho', len(grdROMS.s_rho))
     f1.createDimension('s_w', len(grdROMS.s_w))

     vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
     vnc.long_name = 'Longitude at RHO points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_rho
     
     vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
     vnc.long_name = 'Latitude at RHO points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_rho
     
     vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',),zlib=True)
     vnc.long_name = 'Longitude at U points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_u
     
     vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',),zlib=True)
     vnc.long_name = 'Latitude at U points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_u
     
     vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',),zlib=True)
     vnc.long_name = 'Longitude at V points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_v
     
     vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',),zlib=True)
     vnc.long_name = 'Latitude at V points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_v
     
     vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',),zlib=True)
     vnc.long_name = 'Final bathymetry at RHO points'
     vnc.units = 'meter'
     vnc.field = "bath, scalar"
     vnc[:,:] = grdROMS.depth
     
     vnc = f1.createVariable('s_rho', 'd', ('s_rho',),zlib=True)
     vnc.long_name = "S-coordinate at RHO-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     vnc.standard_name = "ocean_s_coordinate"
     vnc.formula_terms = "s: s_w eta: zeta depth: h a: theta_s b: theta_b depth_c: hc" 
     vnc.field = "s_rho, scalar"
     vnc[:] = grdROMS.s_rho
        
     vnc= f1.createVariable('Cs_rho', 'd', ('s_rho',),zlib=True)
     vnc.long_name = "S-coordinate stretching curves at RHO-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     vnc.field = "s_rho, scalar"
     vnc[:] = grdROMS.Cs_rho
     
     vnc=f1.createVariable('hc','d',zlib=True)
     vnc.long_name = "S-coordinate parameter, critical depth" ;
     vnc.units = "meter"
     vnc[:]=grdROMS.hc
     
     vnc=f1.createVariable('z_r','d',('s_rho','eta_rho','xi_rho',),zlib=True)
     vnc.long_name = "Sigma layer to depth matrix" ;
     vnc.units = "meter"
     vnc[:,:,:]=grdROMS.z_r
    
     vnc=f1.createVariable('Tcline','d',zlib=True)
     vnc.long_name = "S-coordinate surface/bottom layer width" ;
     vnc.units = "meter"
     vnc[:]=grdROMS.Tcline
     
     vnc=f1.createVariable('theta_s','d',zlib=True)
     vnc.long_name = "S-coordinate surface control parameter"
     vnc[:]=grdROMS.theta_s
     
     vnc=f1.createVariable('theta_b','d',zlib=True)
     vnc.long_name = "S-coordinate bottom control parameter"
     vnc[:]=grdROMS.theta_b
     
     vnc=f1.createVariable('angle','d',('eta_rho','xi_rho'),zlib=True)
     vnc.long_name = "angle between xi axis and east"
     vnc.units = "radian" 

     vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho'),zlib=True)
     vnc.long_name = "mask on RHO-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_rho
     
     vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u'),zlib=True)
     vnc.long_name = "mask on U-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_u
     
     vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v'),zlib=True)
     vnc.long_name = "mask on V-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_v
     
     
     v_time = f1.createVariable('bry_time', 'd', ('bry_time',),zlib=True)
     v_time.long_name = 'Days since 1948-01-01 00:00:00'
     v_time.units = 'days'
     v_time.field = 'time, scalar, series'
     v_time.calendar='standard'
     
     v_temp_west=f1.createVariable('temp_west', 'f', ('bry_time', 's_rho', 'eta_rho'),zlib=True)
     v_temp_west.long_name = "Ocean temperature at Western boundary"
     v_temp_west.units = "degrees Celsius"
     v_temp_west.FillValue = grdROMS.fill_value
     v_temp_west.time = "bry_time"
     
     v_temp_east=f1.createVariable('temp_east', 'f', ('bry_time', 's_rho', 'eta_rho'),zlib=True)
     v_temp_east.long_name = "Ocean temperature at Eastern boundary"
     v_temp_east.units = "degrees Celsius"
     v_temp_east.FillValue = grdROMS.fill_value
     v_temp_east.time = "bry_time"
     
     v_temp_south=f1.createVariable('temp_south', 'f', ('bry_time', 's_rho', 'xi_rho'),zlib=True)
     v_temp_south.long_name = "Ocean temperature at Southern boundary"
     v_temp_south.units = "degrees Celsius"
     v_temp_south.FillValue = grdROMS.fill_value
     v_temp_south.time = "bry_time"
     
     v_temp_north=f1.createVariable('temp_north', 'f', ('bry_time', 's_rho', 'xi_rho'),zlib=True)
     v_temp_north.long_name = "Ocean temperature at Northern boundary"
     v_temp_north.units = "degrees Celsius"
     v_temp_north.FillValue = grdROMS.fill_value
     v_temp_north.time = "bry_time"
     
     v_salt_west=f1.createVariable('salt_west', 'f', ('bry_time', 's_rho', 'eta_rho'),zlib=True)
     v_salt_west.long_name = "Ocean salinity at Western boundary"
     v_salt_west.units = "psu"
     v_salt_west.FillValue = grdROMS.fill_value
     v_salt_west.time = "bry_time"
     
     v_salt_east=f1.createVariable('salt_east', 'f', ('bry_time', 's_rho', 'eta_rho'),zlib=True)
     v_salt_east.long_name = "Ocean salinity at Eastern boundary"
     v_salt_east.units = "psu"
     v_salt_east.FillValue = grdROMS.fill_value
     v_salt_east.time = "bry_time"
     
     v_salt_south=f1.createVariable('salt_south', 'f', ('bry_time', 's_rho', 'xi_rho'),zlib=True)
     v_salt_south.long_name = "Ocean salinity at Southern boundary"
     v_salt_south.units = "psu"
     v_salt_south.FillValue = grdROMS.fill_value
     v_salt_south.time = "bry_time"
     
     v_salt_north=f1.createVariable('salt_north', 'f', ('bry_time', 's_rho', 'xi_rho'),zlib=True)
     v_salt_north.long_name = "Ocean salinity at Northern boundary"
     v_salt_north.units = "psu"
     v_salt_north.FillValue = grdROMS.fill_value
     v_salt_north.time = "bry_time"
     
     v_ssh_west=f1.createVariable('zeta_west','d',('bry_time','eta_rho'),zlib=True)
     v_ssh_west.long_name = "Sea surface height (SSH) at Western boundary conditions"
     v_ssh_west.units = "m"
     v_ssh_west.field = "zeta_west, scalar, series"
     v_ssh_west.FillValue = grdROMS.fill_value
     v_ssh_west.time = "bry_time"

     v_ssh_east=f1.createVariable('zeta_east','d',('bry_time','eta_rho'),zlib=True)
     v_ssh_east.long_name = "Sea surface height (SSH) at Eastern boundary conditions"
     v_ssh_east.units = "m"
     v_ssh_east.field = "zeta_east, scalar, series"
     v_ssh_east.FillValue = grdROMS.fill_value
     v_ssh_east.time = "bry_time"
  
     v_ssh_south=f1.createVariable('zeta_south','d',('bry_time','xi_rho'),zlib=True)
     v_ssh_south.long_name = "Sea surface height (SSH) at Southern boundary conditions"
     v_ssh_south.units = "m"
     v_ssh_south.field = "zeta_south, scalar, series"
     v_ssh_south.FillValue = grdROMS.fill_value
     v_ssh_south.time = "bry_time"

     v_ssh_north=f1.createVariable('zeta_north','d',('bry_time','xi_rho'),zlib=True)
     v_ssh_north.long_name = "Sea surface height (SSH) at Northern boundary conditions"
     v_ssh_north.units = "m"
     v_ssh_north.field = "zeta_north, scalar, series"
     v_ssh_north.FillValue = grdROMS.fill_value
     v_ssh_north.time = "bry_time"

     v_u_west=f1.createVariable('u_west', 'f', ('bry_time', 's_rho', 'eta_u'),zlib=True)
     v_u_west.long_name = "U-velocity at Western boundary, scalar, series"
     v_u_west.units = "m/s"
     v_u_west.FillValue = grdROMS.fill_value
     v_u_west.time = "bry_time"
     
     v_u_east=f1.createVariable('u_east', 'f', ('bry_time', 's_rho', 'eta_u'),zlib=True)
     v_u_east.long_name = "U-velocity at Eastern boundary, scalar, series"
     v_u_east.units = "m/s"
     v_u_east.FillValue = grdROMS.fill_value
     v_u_east.time = "bry_time"
     
     v_u_south=f1.createVariable('u_south', 'f', ('bry_time', 's_rho', 'xi_u'),zlib=True)
     v_u_south.long_name = "U-velocity at Southern boundary, scalar, series"
     v_u_south.units = "m/s"
     v_u_south.FillValue = grdROMS.fill_value
     v_u_south.time = "bry_time"
     
     v_u_north=f1.createVariable('u_north', 'f', ('bry_time', 's_rho', 'xi_u'),zlib=True)
     v_u_north.long_name = "U-velocity at Northern boundary, scalar, series"
     v_u_north.units = "m/s"
     v_u_north.FillValue = grdROMS.fill_value
     v_u_north.time = "bry_time"
     
     v_v_west=f1.createVariable('v_west', 'f', ('bry_time', 's_rho', 'eta_v'),zlib=True)
     v_v_west.long_name = "V-velocity at Western boundary, scalar, series"
     v_v_west.units = "m/s"
     v_v_west.FillValue = grdROMS.fill_value
     v_v_west.time = "bry_time"
     
     v_v_east=f1.createVariable('v_east', 'f', ('bry_time', 's_rho', 'eta_v'),zlib=True)
     v_v_east.long_name = "V-velocity at Eastern boundary, scalar, series"
     v_v_east.units = "m/s"
     v_v_east.FillValue = grdROMS.fill_value
     v_v_east.time = "bry_time"
     
     v_v_south=f1.createVariable('v_south', 'f', ('bry_time', 's_rho', 'xi_v'),zlib=True)
     v_v_south.long_name = "V-velocity at Southern boundary, scalar, series"
     v_v_south.units = "m/s"
     v_v_south.FillValue = grdROMS.fill_value
     v_v_south.time = "bry_time"
     
     v_v_north=f1.createVariable('v_north', 'f', ('bry_time', 's_rho', 'xi_v'),zlib=True)
     v_v_north.long_name = "V-velocity at Northern boundary, scalar, series"
     v_v_north.units = "m/s"
     v_v_north.FillValue = grdROMS.fill_value
     v_v_north.time = "bry_time"
     
     v_vbar_west=f1.createVariable('vbar_west', 'f', ('bry_time', 'eta_v'),zlib=True)
     v_vbar_west.long_name = "Barotropic V-velocity at Western boundary, scalar, series"
     v_vbar_west.units = "m/s"
     v_vbar_west.FillValue = grdROMS.fill_value
     v_vbar_west.time = "bry_time"
     
     v_vbar_east=f1.createVariable('vbar_east', 'f', ('bry_time', 'eta_v'),zlib=True)
     v_vbar_east.long_name = "Barotropic V-velocity at Eastern boundary, scalar, series"
     v_vbar_east.units = "m/s"
     v_vbar_east.FillValue = grdROMS.fill_value
     v_vbar_east.time = "bry_time"
     
     v_vbar_south=f1.createVariable('vbar_south', 'f', ('bry_time', 'xi_v'),zlib=True)
     v_vbar_south.long_name = "Barotropic V-velocity at Southern boundary, scalar, series"
     v_vbar_south.units = "m/s"
     v_vbar_south.FillValue = grdROMS.fill_value
     v_vbar_south.time = "bry_time"
     
     v_vbar_north=f1.createVariable('vbar_north', 'f', ('bry_time', 'xi_v'),zlib=True)
     v_vbar_north.long_name = "Barotropic V-velocity at Northern boundary, scalar, series"
     v_vbar_north.units = "m/s"
     v_vbar_north.FillValue = grdROMS.fill_value
     v_vbar_north.time = "bry_time"
     
     v_ubar_west=f1.createVariable('ubar_west', 'f', ('bry_time', 'eta_u'),zlib=True)
     v_ubar_west.long_name = "Barotropic U-velocity at Western boundary, scalar, series"
     v_ubar_west.units = "m/s"
     v_ubar_west.FillValue = grdROMS.fill_value
     v_ubar_west.time = "bry_time"
     
     v_ubar_east=f1.createVariable('ubar_east', 'f', ('bry_time', 'eta_u'),zlib=True)
     v_ubar_east.long_name = "Barotropic U-velocity at Eastern boundary, scalar, series"
     v_ubar_east.units = "m/s"
     v_ubar_east.FillValue = grdROMS.fill_value
     v_ubar_east.time = "bry_time"
     
     v_ubar_south=f1.createVariable('ubar_south', 'f', ('bry_time', 'xi_u'),zlib=True)
     v_ubar_south.long_name = "Barotropic U-velocity at Southern boundary, scalar, series"
     v_ubar_south.units = "m/s"
     v_ubar_south.FillValue = grdROMS.fill_value
     v_ubar_south.time = "bry_time"
     
     v_ubar_north=f1.createVariable('ubar_north', 'f', ('bry_time', 'xi_u'),zlib=True)
     v_ubar_north.long_name = "Barotropic U-velocity at Northern boundary, scalar, series"
     v_ubar_north.units = "m/s"
     v_ubar_north.FillValue = grdROMS.fill_value
     v_ubar_north.time = "bry_time"
     
     f1.close()

