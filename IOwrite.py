
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

def writeClimFile(grdROMS,ntime,outfilename,var,data1=None,data2=None,data3=None,data4=None):
        
     if grdROMS.ioClimInitialized is False:
        grdROMS.ioClimInitialized=True
        if os.path.exists(outfilename):
            os.remove(outfilename)
        
        f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
        f1.title       ="Climatology forcing file (CLIM) used for foring of the ROMS model"
        f1.description ="Created for the %s grid file"%(grdROMS.grdName)
        f1.grd_file    ="Gridfile: %s"%(grdROMS.grdfilename)
        f1.history     ="Created " + time.ctime(time.time())
        f1.source      = "Trond Kristiansen (trond.kristiansen@imr.no)"
        f1.type        ="File in NetCDF4 format created using SODA2ROMS"
        
        """Define dimensions"""   
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
        vnc.long_name = 'Latitude of U-points'
        vnc.units = 'degrees north'
        vnc[:,:] = grdROMS.lat_u
        
        vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',))
        vnc.long_name = 'Longitude of V-points'
        vnc.units = 'degrees east'
        vnc[:,:] = grdROMS.lon_v
        
        vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',))
        vnc.long_name = 'Latitude of V-points'
        vnc.units = 'degrees north'
        vnc[:,:] = grdROMS.lat_v
        
        vnc = f1.createVariable('lat_psi', 'd', ('eta_psi','xi_psi',))
        vnc.long_name = 'Latitude of PSI-points'
        vnc.units = 'degrees north'
        vnc[:,:] = grdROMS.lat_psi
        
        vnc = f1.createVariable('lon_psi', 'd', ('eta_psi','xi_psi',))
        vnc.long_name = 'Longitude of PSI-points'
        vnc.units = 'degrees east'
        vnc[:,:] = grdROMS.lon_psi
        
        vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',))
        vnc.long_name = 'Bathymetry at RHO-points'
        vnc.units = 'meter'
        vnc.field = "bath, scalar"
        vnc[:,:] = grdROMS.depth
        
        vnc = f1.createVariable('f', 'd', ('eta_rho','xi_rho',))
        vnc.long_name = 'Coriolis parameter at RHO-points'
        vnc.units = 'second-1'
        vnc.field = "Coriolis, scalar"
        vnc[:,:] = grdROMS.f
        
        vnc = f1.createVariable('pm', 'd', ('eta_rho','xi_rho',))
        vnc.long_name = 'curvilinear coordinate metric in XI'
        vnc.units = 'meter-1'
        vnc.field = "pm, scalar"
        vnc[:,:] = grdROMS.pm
        
        vnc = f1.createVariable('pn', 'd', ('eta_rho','xi_rho',))
        vnc.long_name = 'curvilinear coordinate metric in ETA'
        vnc.units = 'meter-1'
        vnc.field = "pn, scalar"
        vnc[:,:] = grdROMS.pn
        
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
        vnc[:] = np.flipud(grdROMS.s_rho)
        
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
        vnc[:] = np.flipud(grdROMS.s_w)
        
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
        
        vnc=f1.createVariable('z_w','d',('s_w','eta_rho','xi_rho',))
        vnc.long_name = "Sigma layer to depth matrix" ;
        vnc.units = "meter"
        vnc[:,:,:]=grdROMS.z_w
       
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
        vnc[:,:]=grdROMS.angle
   
        vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho',))
        vnc.long_name = "mask on RHO-points"
        vnc.option_0 = "land" 
        vnc.option_1 = "water"
        vnc.FillValue = 1.0
        vnc[:,:]=grdROMS.mask_rho
        
        vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u',))
        vnc.long_name = "mask at U-points"
        vnc.option_0 = "land" 
        vnc.option_1 = "water"
        vnc.FillValue = 1.0
        vnc[:,:]=grdROMS.mask_u
        
        vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v',))
        vnc.long_name = "mask at V-points"
        vnc.option_0 = "land" 
        vnc.option_1 = "water"
        vnc.FillValue = 1.0
        vnc[:,:]=grdROMS.mask_v
        
        vnc=f1.createVariable('mask_psi','d',('eta_psi', 'xi_psi',))
        vnc.long_name = "mask at PSI-points"
        vnc.option_0 = "land" 
        vnc.option_1 = "water"
        vnc.FillValue = 1.0
        vnc[:,:]=grdROMS.mask_psi
        
        v_time = f1.createVariable('ocean_time', 'd', ('ocean_time',))
        v_time.long_name = 'seconds since 1948-01-01 00:00:00'
        v_time.units = 'seconds since 1948-01-01 00:00:00'
        v_time.field = 'time, scalar, series'
        v_time.calendar='standard'
        
        v_u=f1.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u',))
        v_u.long_name = "u-momentum component"
        v_u.units = "meter second-1"
        v_u.time = "ocean_time"
        v_u.field= "u-velocity, scalar, series"
        v_u.FillValue = grdROMS.fill_value
        
        v_v=f1.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v',))
        v_v.long_name = "v-momentum component"
        v_v.units = "meter second-1"
        v_v.time = "ocean_time"
        v_v.field= "v-velocity, scalar, series"
        v_v.FillValue = grdROMS.fill_value
        
        v_salt=f1.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
        v_salt.long_name = "salinity"
        v_salt.units = "nondimensional"
        v_salt.time = "ocean_time"
        v_salt.field="salinity, scalar, series"
        v_salt.FillValue = grdROMS.fill_value
        
        v_temp=f1.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
        v_temp.long_name = "potential temperature"
        v_temp.units = "Celsius"
        v_temp.time = "ocean_time"
        v_temp.field="temperature, scalar, series"
        v_temp.FillValue = grdROMS.fill_value
        
        v_ssh=f1.createVariable('zeta','f',('ocean_time','eta_rho', 'xi_rho',))
        v_ssh.long_name = "sea level"
        v_ssh.units = "meter"
        v_ssh.time = "ocean_time"
        v_ssh.field="sea level, scalar, series"
        v_ssh.FillValue = grdROMS.fill_value
        
        v_ubar=f1.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u',))
        v_ubar.long_name = "u-2D momentum"
        v_ubar.units = "meter second-1"
        v_ubar.time = "ocean_time"
        v_ubar.field= "u2-D velocity, scalar, series"
        v_ubar.FillValue = grdROMS.fill_value
        
        v_vbar=f1.createVariable('vbar', 'f', ('ocean_time', 'eta_v','xi_v',))
        v_vbar.long_name = "v-2D momentum"
        v_vbar.units = "meter second-1"
        v_vbar.time = "ocean_time"
        v_vbar.field= "v2-D velocity, scalar, series"
        v_vbar.FillValue = grdROMS.fill_value
     
     
     else:
        f1 = Dataset(outfilename, mode='a', format='NETCDF3_CLASSIC')
    
     if var==grdROMS.vars[0]:
          f1.variables['ocean_time'][ntime]   = grdROMS.time * 86400.0
     
     d= num2date(grdROMS.time*86400.0,units=f1.variables['ocean_time'].long_name,calendar=f1.variables['ocean_time'].calendar)
     grdROMS.message=d
     
     if var=='temperature':
        f1.variables['temp'][ntime,:,:,:]  = data1
     if var=='salinity':
        f1.variables['salt'][ntime,:,:,:]  = data1
     if var=='ssh':
        f1.variables['zeta'][ntime,:,:]    = data1
     if var=='vvel':
        
        f1.variables['u'][ntime,:,:,:]     = data1
        f1.variables['v'][ntime,:,:,:]     = data2
     
        f1.variables['ubar'][ntime,:,:]    = data3
        f1.variables['vbar'][ntime,:,:]    = data4
     
     f1.close()


