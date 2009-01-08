
from netCDF4 import Dataset
import numpy as np
import time

def open_output(grdROMS,ntimes):
        
     outfile='test.nc'
     f1 = Dataset(outfile, mode='w', format='NETCDF3_CLASSIC')
     f1.description="This is a climatology file for ROMS"
     f1.history = 'Created ' + time.ctime(time.time())
     f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
     f1.type='NetCDF3 classic created using SODA2ROMS' 
        
     # Define dimensions
     f1.createDimension('xi_rho',  grdROMS.xi_rho)
     f1.createDimension('eta_rho', grdROMS.eta_rho)
     f1.createDimension('xi_u',    grdROMS.xi_u)
     f1.createDimension('eta_u',   grdROMS.eta_u)
     f1.createDimension('xi_v',    grdROMS.xi_v)
     f1.createDimension('eta_v',   grdROMS.eta_v)
     f1.createDimension('xi_psi',    grdROMS.xi_psi)
     f1.createDimension('eta_psi',   grdROMS.eta_psi)
     f1.createDimension('time', ntimes)
     f1.createDimension('s_rho', len(grdROMS.s_rho))
     f1.createDimension('s_w', len(grdROMS.s_w))
     
     v = f1.createVariable('clim_time', 'd', ('time',))
     v.long_name = 'Time since 1948/01/01/00:00'
     v.units = 'days'
     v.field = '"time, scalar, series"'
     v[:] = grdROMS.time
     
     v = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',))
     v.long_name = 'Longitude at RHO points'
     v.units = 'degrees east'
     v[:,:] = grdROMS.lon_rho
     
     
     v = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',))
     v.long_name = 'Latitude at RHO points'
     v.units = 'degrees north'
     v[:,:] = grdROMS.lat_rho
     
     v = f1.createVariable('h', 'd', ('eta_rho','xi_rho',))
     v.long_name = 'Final bathymetry at RHO points'
     v.units = 'meter'
     v.field = "bath, scalar"
     v[:,:] = grdROMS.depth
     
     v = f1.createVariable('s_rho', 'd', ('s_rho',))
     v0 = grdROMS.s_rho
     v.long_name = "S-coordinate at RHO-points"
     v.valid_min = -1.
     v.valid_max = 0.
     v.standard_name = "ocean_s_coordinate"
     v.formula_terms = "s: s_w eta: zeta depth: h a: theta_s b: theta_b depth_c: hc" 
     v.field = "s_rho, scalar"
     v[:] = v0[:]
        
     v  = f1.createVariable('Cs_rho', 'd', ('s_rho',))
     v0 = grdROMS.Cs_rho
     v.long_name = "S-coordinate stretching curves at RHO-points"
     v.valid_min = -1.
     v.valid_max = 0.
     v.field = "s_rho, scalar"
     v[:] = v0[:]
     
     v=f1.createVariable('hc','d')
     v.long_name = "S-coordinate parameter, critical depth" ;
     v.units = "meter"
     v[:]=grdROMS.hc
     
     v=f1.createVariable('z_r','d',('s_rho','eta_rho','xi_rho',))
     v.long_name = "Sigma layer to depth matrix" ;
     v.units = "meter"
     v[:,:,:]=grdROMS.z_r
    
     v=f1.createVariable('Tcline','d')
     v.long_name = "S-coordinate surface/bottom layer width" ;
     v.units = "meter"
     v[:]=grdROMS.Tcline
     
     v=f1.createVariable('theta_s','d')
     v.long_name = "S-coordinate surface control parameter"
     v[:]=grdROMS.theta_s
     
     v=f1.createVariable('theta_b','d')
     v.long_name = "S-coordinate bottom control parameter"
     v[:]=grdROMS.theta_b
     
     v=f1.createVariable('angle','d',('eta_rho','xi_rho'))
     v.long_name = "angle between xi axis and east"
     v.units = "radian" 

     v=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho'))
     v.long_name = "mask on RHO-points"
     v.option_0 = "land" 
     v.option_1 = "water"
     v.FillValue = 1.0
     v[:,:]=grdROMS.mask_rho
     
     v=f1.createVariable('mask_u','d',('eta_u', 'xi_u'))
     v.long_name = "mask on U-points"
     v.option_0 = "land" 
     v.option_1 = "water"
     v.FillValue = 1.0
     v[:,:]=grdROMS.mask_u
     
     v=f1.createVariable('mask_v','d',('eta_v', 'xi_v'))
     v.long_name = "mask on V-points"
     v.option_0 = "land" 
     v.option_1 = "water"
     v.FillValue = 1.0
     v[:,:]=grdROMS.mask_v
     
     v=f1.createVariable('temp', 'f', ('time', 's_rho', 'eta_rho', 'xi_rho'))
     v.long_name = "Ocean temperature"
     v.units = "degrees Celsius"
     v.FillValue = grdROMS.fill_value
     v[:,:,:,:]=grdROMS.t2
     
     v=f1.createVariable('salt', 'f', ('time', 's_rho', 'eta_rho', 'xi_rho'))
     v.long_name = "Ocean salinity"
     v.units = "psu"
     v.FillValue = grdROMS.fill_value
     v[:,:,:,:]=grdROMS.s2
     
     f1.close()


