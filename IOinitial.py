import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 3, 16)
__modified__ = datetime(2009, 11, 11)
__version__  = "1.0"
__status__   = "Development"


def help ():
     """
     This function generates an INIT file from scratch. Varibales include:
     salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is ocean_time which is days
     since 1948/1/1.
     
     Edited by Trond Kristiansen, 16.3.2009, 11.11.2009, 20.11.2009
     """

def createInitFile(grdROMS,ntime,outfilename,var,data1=None,data2=None,data3=None,data4=None):
     
     """
     Create initial file for use with ROMS. This is the same as extracting time 0 from
     the climatology file.
     """
     if grdROMS.ioInitInitialized is False:
          grdROMS.ioInitInitialized=True
          if os.path.exists(outfilename):
               os.remove(outfilename)
       
          f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
          f1.title       ="Initial forcing file (INIT) used for foring of the ROMS model"
          f1.description ="Created for the %s grid file"%(grdROMS.grdName)
          f1.grd_file    ="Gridfile: %s"%(grdROMS.grdfilename)
          f1.history     ="Created " + time.ctime(time.time())
          f1.source      ="Trond Kristiansen (trond.kristiansen@imr.no)"
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
          vnc.long_name = 'Longitude at RHO points'
          vnc.units = 'degrees east'
          vnc[:,:] = grdROMS.lon_rho
          
          vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',))
          vnc.long_name = 'Latitude at RHO points'
          vnc.units = 'degrees north'
          vnc[:,:] = grdROMS.lat_rho
          
          vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',))
          vnc.long_name = 'Longitude at U points'
          vnc.units = 'degrees east'
          vnc[:,:] = grdROMS.lon_u
          
          vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',))
          vnc.long_name = 'Latitude at U points'
          vnc.units = 'degrees north'
          vnc[:,:] = grdROMS.lat_u
          
          vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',))
          vnc.long_name = 'Longitude at V points'
          vnc.units = 'degrees east'
          vnc[:,:] = grdROMS.lon_v
          
          vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',))
          vnc.long_name = 'Latitude at V points'
          vnc.units = 'degrees north'
          vnc[:,:] = grdROMS.lat_v
          
          vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',))
          vnc.long_name = 'Final bathymetry at RHO points'
          vnc.units = 'meter'
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
        
          vnc= f1.createVariable('Cs_rho', 'd', ('s_rho',))
          vnc.long_name = "S-coordinate stretching curves at RHO-points"
          vnc.valid_min = -1.
          vnc.valid_max = 0.
          vnc.field = "s_rho, scalar"
          vnc[:] = grdROMS.Cs_rho
          
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
          v_time.units = 'seconds'
          v_time.field = 'time, scalar, series'
          v_time.calendar='standard'
          
          v_temp=f1.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
          v_temp.long_name = "potential temperature"
          v_temp.units = "Celsius"
          v_temp.time = "ocean_time"
          v_temp.FillValue = grdROMS.fill_value
          
          v_salt=f1.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
          v_salt.long_name = "salinity"
          v_salt.units = "PSU"
          v_salt.time = "ocean_time"
          v_salt.FillValue = grdROMS.fill_value
          
          v_ssh=f1.createVariable('zeta','f',('ocean_time','eta_rho', 'xi_rho',))
          v_ssh.long_name = "sea level"
          v_ssh.units = "meter"
          v_ssh.time = "ocean_time"
          v_ssh.FillValue = grdROMS.fill_value
          
          v_u=f1.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u',))
          v_u.long_name = "U-velocity, scalar, series"
          v_u.units = "meter second-1"
          v_u.time = "ocean_time"
          v_u.FillValue = grdROMS.fill_value
          
          v_v=f1.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v',))
          v_v.long_name = "V-velocity, scalar, series"
          v_v.units = "meter second-1"
          v_v.time = "ocean_time"
          v_v.FillValue = grdROMS.fill_value
          
          v_vbar=f1.createVariable('vbar', 'f', ('ocean_time', 'eta_v','xi_v',))
          v_vbar.long_name = "Barotropic V-velocity, scalar, series"
          v_vbar.units = "meter second-1"
          v_vbar.time = "ocean_time"
          v_vbar.FillValue = grdROMS.fill_value
        
          v_ubar=f1.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u',))
          v_ubar.long_name = "Barotropic U-velocity, scalar, series"
          v_ubar.units = "meter second-1"
          v_ubar.time = "ocean_time"
          v_ubar.FillValue = grdROMS.fill_value
              
          d= num2date(grdROMS.time,units=v_time.long_name,calendar=v_time.calendar)
          print '\n'
          print '========================================================================='
          print 'Created INIT file'
          print 'set initTime in grd.py for time-index to print (current=%s)'%(grdROMS.initTime)
          print 'The time stamp for ROMS .in file saved to initial file is=> %s'%(d)
          print 'DSTART   = %s'%(grdROMS.time)
          print 'TIME_REF = %s'%(v_time.long_name)
          print '========================================================================='
          print '\n'
     
     else:
        f1 = Dataset(outfilename, mode='a', format='NETCDF4', zlib=True)
        
     ntime=0
     f1.variables['ocean_time'][ntime]   = grdROMS.time * 86400.0
       
     if var=='temperature':
        f1.variables['temp'][ntime,:,:,:]  = data1
     if var=='salinity':
        f1.variables['salt'][ntime,:,:,:]  = data1
     if var=='ssh':
        f1.variables['zeta'][ntime,:,:]    = data1
     if var in ['uvel','vvel','ubar','vbar']:
        f1.variables['u'][ntime,:,:,:]     = data1
        f1.variables['v'][ntime,:,:,:]     = data2
        f1.variables['ubar'][ntime,:,:]    = data3
        f1.variables['vbar'][ntime,:,:]    = data4
     
     f1.close()


