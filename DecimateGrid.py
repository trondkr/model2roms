
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

def createGrid(grdROMS,outfilename,decimate):
        
        
    if os.path.exists(outfilename):
        os.remove(outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a grid file for ROMS"
    f1.history = 'Created (decimated) from IMR Nordic 4km grid' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 classic created using SODA2ROMS'
    
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_rho),int(grdROMS.eta_rho))
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_u),int(grdROMS.eta_u))
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_v),int(grdROMS.eta_v))
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_psi),int(grdROMS.eta_psi))
     
    print 'New dimensions will be : %ix%i'%(int(grdROMS.xi_rho/2.0),int(grdROMS.eta_rho/2.0))
    print 'New dimensions will be : %ix%i'%(int(grdROMS.xi_u/2.0),int(grdROMS.eta_u/2.0))
    print 'New dimensions will be : %ix%i'%(int(grdROMS.xi_v/2.0),int(grdROMS.eta_v/2.0))
    print 'New dimensions will be : %ix%i'%(int(grdROMS.xi_psi/2.0),int(grdROMS.eta_psi/2.0))
    
    # Define dimensions
    f1.createDimension('xi_rho',  int(grdROMS.xi_rho/decimate))
    f1.createDimension('eta_rho', int(grdROMS.eta_rho/decimate))
    f1.createDimension('xi_u',    int(grdROMS.xi_u/decimate))
    f1.createDimension('eta_u',   int(grdROMS.eta_u/decimate))
    f1.createDimension('xi_v',    int(grdROMS.xi_v/decimate))
    f1.createDimension('eta_v',   int(grdROMS.eta_v/decimate))
    f1.createDimension('xi_psi',    int(grdROMS.xi_psi/decimate))
    f1.createDimension('eta_psi',   int(grdROMS.eta_psi/decimate))
    f1.createDimension('s_rho', int(len(grdROMS.s_rho)/decimate))
    f1.createDimension('s_w', int(len(grdROMS.s_w)/decimate))

    vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Longitude at RHO points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_rho[::decimate,::decimate]
    
    vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Latitude at RHO points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_rho[::decimate,::decimate]
    
    vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',),zlib=True)
    vnc.long_name = 'Longitude at U points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_u[::decimate,0:-1:decimate]
    
    vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',),zlib=True)
    vnc.long_name = 'Latitude at U points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_u[::decimate,0:-1:decimate]
    
    print grdROMS.lon_v[::decimate,0:-1:decimate].shape
    vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',),zlib=True)
    vnc.long_name = 'Longitude at V points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_v[0:-1:decimate,::decimate]
    
    vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',),zlib=True)
    vnc.long_name = 'Latitude at V points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_v[0:-1:decimate,::decimate]
    
    vnc = f1.createVariable('lat_psi', 'd', ('eta_psi','xi_psi',),zlib=True)
    vnc.long_name = 'Latitude at PSI points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_psi[0:-1:decimate,0:-1:decimate]
    
    vnc = f1.createVariable('lon_psi', 'd', ('eta_psi','xi_psi',),zlib=True)
    vnc.long_name = 'Longitude at PSI points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_psi[0:-1:decimate,0:-1:decimate]
    
    vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Final bathymetry at RHO points'
    vnc.units = 'meter'
    vnc.field = "bath, scalar"
    vnc[:,:] = grdROMS.depth[::decimate,::decimate]
    
    vnc = f1.createVariable('f', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Coriolis parameter at RHO points'
    vnc.units = 'second-1'
    vnc.field = "Coriolis, scalar"
    vnc[:,:] = grdROMS.f[::decimate,::decimate]
    
    vnc = f1.createVariable('pm', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'curvilinear coordinate metric in XI'
    vnc.units = 'meter-1'
    vnc.field = "pm, scalar"
    vnc[:,:] = grdROMS.pm[::decimate,::decimate]
    
    vnc = f1.createVariable('pn', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'curvilinear coordinate metric in ETA'
    vnc.units = 'meter-1'
    vnc.field = "pn, scalar"
    vnc[:,:] = grdROMS.pn[::decimate,::decimate]
    
    
    vnc=f1.createVariable('angle','d',('eta_rho','xi_rho'),zlib=True)
    vnc.long_name = "angle between xi axis and east"
    vnc.units = "radian" 
    vnc[:,:]=grdROMS.angle[::decimate,::decimate]
    
    vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho'),zlib=True)
    vnc.long_name = "mask on RHO-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_rho[::decimate,::decimate]

    
    vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u'),zlib=True)
    vnc.long_name = "mask on U-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_u[::decimate,0:-1:decimate]
    
    vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v'),zlib=True)
    vnc.long_name = "mask on V-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_v[0:-1:decimate,::decimate]
    
    vnc=f1.createVariable('mask_psi','d',('eta_psi', 'xi_psi'),zlib=True)
    vnc.long_name = "mask on PSI-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_psi[0:-1:decimate,0:-1:decimate]
    
    f1.close()


