import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 3,2)
__modified__ = datetime(2009, 11,18)
__version__  = "0.9.1"
__status__   = "Development"


def help ():
     """
     This function generates a BRY file from scratch. The variables are all created
     for East, West, North, and South. Varibales include:
     salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is ocean_time which is days
     since 1948/1/1.
     
     This file is netcdf CF compliant and follows the setup for vairable names and units given in the ROMS source
     file: Data/ROMS/CDL/bry_unlimit.cdl
     
     (also see: https://www.myroms.org/forum/viewtopic.php?f=30&t=1450&p=5209&hilit=cf+compliant#p5209)

     This function is called from clim2bry.py.
     """

def createBryFile(grdROMS,outfilename):
     
     if os.path.exists(outfilename):
         os.remove(outfilename)
     print 'Creating initial Boundary (BRY) file %s'%(outfilename)
     
     f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
     f1.title="Boundary forcing file (BRY) used for foring of the ROMS model"
     f1.description="Created for the %s grid file"%(grdROMS.grdName)
     f1.grdFile="%s"%(grdROMS.grdfilename)
     f1.history = 'Created ' + time.ctime(time.time())
     f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
     f1.type='File in NetCDF4 format created using SODA2ROMS'
     
     
     """ Define dimensions """
     f1.createDimension('xi_rho',  grdROMS.xi_rho)
     f1.createDimension('eta_rho', grdROMS.eta_rho)
     f1.createDimension('xi_u',    grdROMS.xi_u)
     f1.createDimension('eta_u',   grdROMS.eta_u)
     f1.createDimension('xi_v',    grdROMS.xi_v)
     f1.createDimension('eta_v',   grdROMS.eta_v)
     f1.createDimension('xi_psi',    grdROMS.xi_psi)
     f1.createDimension('eta_psi',   grdROMS.eta_psi)
     f1.createDimension('ocean_time', None)
     f1.createDimension('s_rho', len(grdROMS.s_rho))
     f1.createDimension('s_w', len(grdROMS.s_w))

     vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',))
     vnc.long_name = 'Longitude of RHO-points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_rho
     
     vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',))
     vnc.long_name = 'Latitude of RHO-points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_rho
     
     vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',))
     vnc.long_name = 'Longitude of U-points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_u
     
     vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',))
     vnc.long_name = 'Latitude af U-points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_u
     
     vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',))
     vnc.long_name = 'Longitude of V-points'
     vnc.units = 'degrees east'
     vnc[:,:] = grdROMS.lon_v
     
     vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',))
     vnc.long_name = 'Latitude af V-points'
     vnc.units = 'degrees north'
     vnc[:,:] = grdROMS.lat_v
     
     vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',))
     vnc.long_name = "Bathymetry at RHO-points"
     vnc.units = "meter"
     vnc.coordinates ="lat_rho lon_rho"
     vnc.field = "bath, scalar"
     vnc[:,:] = grdROMS.depth
     
     vnc = f1.createVariable('s_rho', 'd', ('s_rho',))
     vnc.long_name = "S-coordinate at RHO-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     if grdROMS.vstretching==2:
          vnc.standard_name = "ocean_s_coordinate_g2" 
          vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
     if grdROMS.vstretching==1:
          vnc.standard_name = "ocean_s_coordinate_g1" 
          vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" 
     vnc.field = "s_rho, scalar"
     vnc[:] = grdROMS.s_rho
        
     vnc = f1.createVariable('s_w', 'd', ('s_w',))
     vnc.long_name = "S-coordinate at W-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     if grdROMS.vstretching==2:
          vnc.standard_name = "ocean_s_coordinate_g2" 
          vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
     if grdROMS.vstretching==1:
          vnc.standard_name = "ocean_s_coordinate_g1" 
          vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc" 
     vnc.field = "s_w, scalar"
     vnc[:] = grdROMS.s_w
     
     vnc= f1.createVariable('Cs_r', 'd', ('s_rho',))
     vnc.long_name = "S-coordinate stretching curves at RHO-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     vnc.field = "s_rho, scalar"
     vnc[:] = np.flipud(grdROMS.Cs_rho)
        
     vnc= f1.createVariable('Cs_w', 'd', ('s_w',))
     vnc.long_name = "S-coordinate stretching curves at W-points"
     vnc.valid_min = -1.
     vnc.valid_max = 0.
     vnc.field = "s_w, scalar"
     vnc[:] = np.flipud(grdROMS.Cs_w)
     
     vnc=f1.createVariable('hc','d')
     vnc.long_name = "S-coordinate parameter, critical depth" ;
     vnc.units = "meter"
     vnc[:]=grdROMS.hc
     
     vnc=f1.createVariable('z_r','d',('s_rho','eta_rho','xi_rho',))
     vnc.long_name = "Sigma layer to depth matrix" ;
     vnc.units = "meter"
     vnc[:,:,:]=grdROMS.z_r
    
     vnc=f1.createVariable('Tcline','d')
     vnc.long_name = "S-coordinate surface/bottom layer width" ;
     vnc.units = "meter"
     vnc[:]=grdROMS.Tcline
     
     vnc=f1.createVariable('theta_s','d')
     vnc.long_name = "S-coordinate surface control parameter"
     vnc[:]=grdROMS.theta_s
     
     vnc=f1.createVariable('theta_b','d')
     vnc.long_name = "S-coordinate bottom control parameter"
     vnc[:]=grdROMS.theta_b
     
     vnc=f1.createVariable('angle','d',('eta_rho','xi_rho',))
     vnc.long_name = "angle between xi axis and east"
     vnc.units = "radian" 

     vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho',))
     vnc.long_name = "mask on RHO-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_rho
     
     vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u',))
     vnc.long_name = "mask on U-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_u
     
     vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v',))
     vnc.long_name = "mask on V-points"
     vnc.option_0 = "land" 
     vnc.option_1 = "water"
     vnc.FillValue = 1.0
     vnc[:,:]=grdROMS.mask_v
     
     v_time = f1.createVariable('ocean_time', 'd', ('ocean_time',))
     v_time.long_name = 'seconds since 1948-01-01 00:00:00'
     v_time.units = 'seconds since 1948-01-01 00:00:00'
     v_time.field = 'time, scalar, series'
     v_time.calendar='standard'
     
     v_temp_west=f1.createVariable('temp_west', 'f', ('ocean_time', 's_rho', 'eta_rho',))
     v_temp_west.long_name = "potential temperature western boundary condition"
     v_temp_west.units = "Celsius"
     v_temp_west.field = "temp_west, scalar, series"
     v_temp_west.FillValue = grdROMS.fill_value
     v_temp_west.time = "ocean_time"
     
     v_temp_east=f1.createVariable('temp_east', 'f', ('ocean_time', 's_rho', 'eta_rho',))
     v_temp_east.long_name = "potential temperature eastern boundary condition"
     v_temp_east.units = "Celsius"
     v_temp_east.field = "temp_east, scalar, series"
     v_temp_east.FillValue = grdROMS.fill_value
     v_temp_east.time = "ocean_time"
     
     v_temp_south=f1.createVariable('temp_south', 'f', ('ocean_time', 's_rho', 'xi_rho',))
     v_temp_south.long_name = "potential temperature southern boundary condition"
     v_temp_south.units = "Celsius"
     v_temp_south.field = "temp_south, scalar, series"
     v_temp_south.FillValue = grdROMS.fill_value
     v_temp_south.time = "ocean_time"
     
     v_temp_north=f1.createVariable('temp_north', 'f', ('ocean_time', 's_rho', 'xi_rho',))
     v_temp_north.long_name = "potential temperature northern boundary condition"
     v_temp_north.units = "Celsius"
     v_temp_north.field = "temp_north, scalar, series"
     v_temp_north.FillValue = grdROMS.fill_value
     v_temp_north.time = "ocean_time"
     
     v_salt_west=f1.createVariable('salt_west', 'f', ('ocean_time', 's_rho', 'eta_rho',))
     v_salt_west.long_name = "salinity western boundary condition"
     v_salt_west.units = "nondimensional"
     v_salt_west.field = "salt_west, scalar, series"
     v_salt_west.FillValue = grdROMS.fill_value
     v_salt_west.time = "ocean_time"
     
     v_salt_east=f1.createVariable('salt_east', 'f', ('ocean_time', 's_rho', 'eta_rho',))
     v_salt_east.long_name = "salinity eastern boundary condition"
     v_salt_east.units = "nondimensional"
     v_salt_east.field = "salt_east, scalar, series"
     v_salt_east.FillValue = grdROMS.fill_value
     v_salt_east.time = "ocean_time"
     
     v_salt_south=f1.createVariable('salt_south', 'f', ('ocean_time', 's_rho', 'xi_rho',))
     v_salt_south.long_name = "salinity southern boundary condition"
     v_salt_south.units = "nondimensional"
     v_salt_south.field = "salt_south, scalar, series"
     v_salt_south.FillValue = grdROMS.fill_value
     v_salt_south.time = "ocean_time"
     
     v_salt_north=f1.createVariable('salt_north', 'f', ('ocean_time', 's_rho', 'xi_rho',))
     v_salt_north.long_name = "salinity northern boundary condition"
     v_salt_north.units = "nondimensional"
     v_salt_north.field = "salt_north, scalar, series"
     v_salt_north.FillValue = grdROMS.fill_value
     v_salt_north.time = "ocean_time"
     
     v_ssh_west=f1.createVariable('zeta_west','f',('ocean_time','eta_rho',))
     v_ssh_west.long_name = "free-surface western boundary condition"
     v_ssh_west.units = "meter"
     v_ssh_west.field = "zeta_west, scalar, series"
     v_ssh_west.FillValue = grdROMS.fill_value
     v_ssh_west.time = "ocean_time"

     v_ssh_east=f1.createVariable('zeta_east','f',('ocean_time','eta_rho',))
     v_ssh_east.long_name = "free-surface eastern boundary condition"
     v_ssh_east.units = "meter"
     v_ssh_east.field = "zeta_east, scalar, series"
     v_ssh_east.FillValue = grdROMS.fill_value
     v_ssh_east.time = "ocean_time"
  
     v_ssh_south=f1.createVariable('zeta_south','f',('ocean_time','xi_rho',))
     v_ssh_south.long_name = "free-surface southern boundary condition"
     v_ssh_south.units = "meter"
     v_ssh_south.field = "zeta_south, scalar, series"
     v_ssh_south.FillValue = grdROMS.fill_value
     v_ssh_south.time = "ocean_time"

     v_ssh_north=f1.createVariable('zeta_north','f',('ocean_time','xi_rho',))
     v_ssh_north.long_name = "free-surface northern boundary condition"
     v_ssh_north.units = "meter"
     v_ssh_north.field = "zeta_north, scalar, series"
     v_ssh_north.FillValue = grdROMS.fill_value
     v_ssh_north.time = "ocean_time"

     v_u_west=f1.createVariable('u_west', 'f', ('ocean_time', 's_rho', 'eta_u',))
     v_u_west.long_name = "3D u-momentum western boundary condition"
     v_u_west.units = "meter second-1"
     v_u_west.field = "u_west, scalar, series"
     v_u_west.FillValue = grdROMS.fill_value
     v_u_west.time = "ocean_time"
     
     v_u_east=f1.createVariable('u_east', 'f', ('ocean_time', 's_rho', 'eta_u',))
     v_u_east.long_name = "3D u-momentum eastern boundary condition"
     v_u_east.units = "meter second-1"
     v_u_east.field = "u_east, scalar, series"
     v_u_east.FillValue = grdROMS.fill_value
     v_u_east.time = "ocean_time"
     
     v_u_south=f1.createVariable('u_south', 'f', ('ocean_time', 's_rho', 'xi_u',))
     v_u_south.long_name = "3D u-momentum southern boundary condition"
     v_u_south.units = "meter second-1"
     v_u_south.field = "u_south, scalar, series"
     v_u_south.FillValue = grdROMS.fill_value
     v_u_south.time = "ocean_time"
     
     v_u_north=f1.createVariable('u_north', 'f', ('ocean_time', 's_rho', 'xi_u',))
     v_u_north.long_name = "3D u-momentum northern boundary condition"
     v_u_north.units = "meter second-1"
     v_u_north.field = "u_north, scalar, series"
     v_u_north.FillValue = grdROMS.fill_value
     v_u_north.time = "ocean_time"
     
     v_v_west=f1.createVariable('v_west', 'f', ('ocean_time', 's_rho', 'eta_v',))
     v_v_west.long_name = "3D v-momentum western boundary condition"
     v_v_west.units = "meter second-1"
     v_v_west.field = "v_west, scalar, series"
     v_v_west.FillValue = grdROMS.fill_value
     v_v_west.time = "ocean_time"
     
     v_v_east=f1.createVariable('v_east', 'f', ('ocean_time', 's_rho', 'eta_v',))
     v_v_east.long_name = "3D v-momentum eastern boundary condition"
     v_v_east.units = "meter second-1"
     v_v_east.field = "v_east, scalar, series"
     v_v_east.FillValue = grdROMS.fill_value
     v_v_east.time = "ocean_time"
     
     v_v_south=f1.createVariable('v_south', 'f', ('ocean_time', 's_rho', 'xi_v',))
     v_v_south.long_name = "3D v-momentum southern boundary condition"
     v_v_south.units = "meter second-1"
     v_v_south.field = "v_south, scalar, series"
     v_v_south.FillValue = grdROMS.fill_value
     v_v_south.time = "ocean_time"
     
     v_v_north=f1.createVariable('v_north', 'f', ('ocean_time', 's_rho', 'xi_v',))
     v_v_north.long_name = "3D v-momentum northern boundary condition"
     v_v_north.units = "meter second-1"
     v_v_north.field = "v_north, scalar, series"
     v_v_north.FillValue = grdROMS.fill_value
     v_v_north.time = "ocean_time"
     
     v_vbar_west=f1.createVariable('vbar_west', 'f', ('ocean_time', 'eta_v',))
     v_vbar_west.long_name = "2D v-momentum western boundary condition"
     v_vbar_west.units = "meter second-1"
     v_vbar_west.field = "vbar_west, scalar, series"
     v_vbar_west.FillValue = grdROMS.fill_value
     v_vbar_west.time = "ocean_time"
     
     v_vbar_east=f1.createVariable('vbar_east', 'f', ('ocean_time', 'eta_v',))
     v_vbar_east.long_name = "2D v-momentum eastern boundary condition"
     v_vbar_east.units = "meter second-1"
     v_vbar_east.field = "vbar_east, scalar, series"
     v_vbar_east.FillValue = grdROMS.fill_value
     v_vbar_east.time = "ocean_time"
     
     v_vbar_south=f1.createVariable('vbar_south', 'f', ('ocean_time', 'xi_v',))
     v_vbar_south.long_name = "2D v-momentum southern boundary condition"
     v_vbar_south.units = "meter second-1"
     v_vbar_south.field = "vbar_south, scalar, series"
     v_vbar_south.FillValue = grdROMS.fill_value
     v_vbar_south.time = "ocean_time"
     
     v_vbar_north=f1.createVariable('vbar_north', 'f', ('ocean_time', 'xi_v',))
     v_vbar_north.long_name = "2D v-momentum northern boundary condition"
     v_vbar_north.units = "meter second-1"
     v_vbar_north.field = "vbar_north, scalar, series"
     v_vbar_north.FillValue = grdROMS.fill_value
     v_vbar_north.time = "ocean_time"
     
     v_ubar_west=f1.createVariable('ubar_west', 'f', ('ocean_time', 'eta_u',))
     v_ubar_west.long_name = "2D u-momentum western boundary condition"
     v_ubar_west.units = "meter second-1"
     v_ubar_west.field = "ubar_west, scalar, series"
     v_ubar_west.FillValue = grdROMS.fill_value
     v_ubar_west.time = "ocean_time"
     
     v_ubar_east=f1.createVariable('ubar_east', 'f', ('ocean_time', 'eta_u',))
     v_ubar_east.long_name = "2D u-momentum eastern boundary condition"
     v_ubar_east.units = "meter second-1"
     v_ubar_east.field = "ubar_east, scalar, series"
     v_ubar_east.FillValue = grdROMS.fill_value
     v_ubar_east.time = "ocean_time"
     
     v_ubar_south=f1.createVariable('ubar_south', 'f', ('ocean_time', 'xi_u',))
     v_ubar_south.long_name = "2D u-momentum southern boundary condition"
     v_ubar_south.units = "meter second-1"
     v_ubar_south.field = "ubar_south, scalar, series"
     v_ubar_south.FillValue = grdROMS.fill_value
     v_ubar_south.time = "ocean_time"
     
     v_ubar_north=f1.createVariable('ubar_north', 'f', ('ocean_time', 'xi_u',))
     v_ubar_north.long_name = "2D u-momentum northern boundary condition"
     v_ubar_north.units = "meter second-1"
     v_ubar_north.field = "ubar_north, scalar, series"
     v_ubar_north.FillValue = grdROMS.fill_value
     v_ubar_north.time = "ocean_time"
     
     f1.close()

