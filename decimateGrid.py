
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os
import shutil
import calculateGRDMetrics

def createGrid(grdROMS,infilename,outfilename,decimate):
        
    shutil.copy2(infilename, outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a grid file for ROMS - KINO project"
    f1.history = 'Created (decimated) from IMR 800M grid' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 classic created using MODEL2ROMS - https://github.com/trondkr/model2roms'
    
    if not int(grdROMS.xi_rho) % 2 == 0: 
        deltaXI=1
    else:
        deltaXI=0
    if not int(grdROMS.eta_rho) % 2 == 0: 
        deltaETA=1
        startindex=0;endindex=-2
    else:
        deltaETA=0
        startindex=0;endindex=-1

    print('Old dimensions were  (XI_RHO, ETA_RHO)   : %ix%i'%(int(grdROMS.xi_rho)-deltaXI,int(grdROMS.eta_rho)-deltaETA))
    print('Old dimensions were  (XI_U, ETA_U)       : %ix%i'%(int(grdROMS.xi_u)-deltaXI,int(grdROMS.eta_u)-deltaETA))
    print('Old dimensions were  (XI_V, ETA_V)       : %ix%i'%(int(grdROMS.xi_v)-deltaXI,int(grdROMS.eta_v)-deltaETA))
    print('Old dimensions were  (XI_PSI, ETA_PSI)   : %ix%i'%(int(grdROMS.xi_psi)-deltaXI,int(grdROMS.eta_psi)-deltaETA))
     
    print('New dimensions will be (XI_RHO, ETA_RHO) : %ix%i'%(((int(grdROMS.xi_rho)-deltaXI)/decimate),((int(grdROMS.eta_rho)-deltaETA)/decimate)))
    print('New dimensions will be (XI_U, ETA_U)     : %ix%i'%(((int(grdROMS.xi_u)-deltaXI)/decimate),((int(grdROMS.eta_u)-deltaETA)/decimate)))
    print('New dimensions will be (XI_V, ETA_V)     : %ix%i'%(((int(grdROMS.xi_v)-deltaXI)/decimate),((int(grdROMS.eta_v)-deltaETA)/decimate)))
    print('New dimensions will be (XI_PSI, ETA_PSI) : %ix%i'%(((int(grdROMS.xi_psi)-deltaXI)/decimate),((int(grdROMS.eta_psi)-deltaETA)/decimate)))
    
    # Define dimensions
    xi_rho=(int(grdROMS.xi_rho)-deltaXI)/decimate
    xi_vert=((int(grdROMS.xi_rho)-deltaXI)/decimate) + 1
    eta_rho=(int(grdROMS.eta_rho)-deltaETA)/decimate
    eta_vert=((int(grdROMS.eta_rho)-deltaETA)/decimate) + 1 
    xi_u=(int(grdROMS.xi_u)-deltaXI)/decimate
    eta_u=(int(grdROMS.eta_u)-deltaETA)/decimate
    xi_v=(int(grdROMS.xi_v)-deltaXI)/decimate
    eta_v=(int(grdROMS.eta_v)-deltaETA)/decimate

    xi_psi=(int(grdROMS.xi_psi)-deltaXI)/decimate
    eta_psi=(int(grdROMS.eta_psi)-deltaETA)/decimate
    s_rho=int(len(grdROMS.s_rho))
    s_w=int(len(grdROMS.s_w))

    f1.createDimension('xi_rho',  xi_rho)
    f1.createDimension('eta_rho', eta_rho)
    f1.createDimension('xi_u',    xi_u)
    f1.createDimension('eta_u',   eta_u)
    f1.createDimension('xi_v',    xi_v)
    f1.createDimension('eta_v',   eta_v)
    f1.createDimension('xi_psi',  xi_psi)
    f1.createDimension('eta_psi', eta_psi)
    f1.createDimension('xi_vert', xi_vert)
    f1.createDimension('eta_vert', eta_vert)
    f1.createDimension('s_rho',   s_rho)
    f1.createDimension('s_w',     s_w)
    f1.createDimension('bath',     None)

    vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Longitude at RHO points'
    vnc.units = 'degree_east'
    vnc[:,:] = grdROMS.lon_rho[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Latitude at RHO points'
    vnc.units = 'degree_north'
    vnc[:,:] = grdROMS.lat_rho[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',),zlib=False)
    vnc.long_name = 'Longitude at U points'
    vnc.units = 'degree_east'
    vnc[:,:] = grdROMS.lon_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',),zlib=False)
    vnc.long_name = 'Latitude at U points'
    vnc.units = 'degree_north'
    vnc[:,:] = grdROMS.lat_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',),zlib=False)
    vnc.long_name = 'Longitude at V points'
    vnc.units = 'degree_east'
    vnc[:,:] = grdROMS.lon_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',),zlib=False)
    vnc.long_name = 'Latitude at V points'
    vnc.units = 'degree_north'
    vnc[:,:] = grdROMS.lat_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_psi', 'd', ('eta_psi','xi_psi',),zlib=False)
    vnc.long_name = 'Latitude at PSI points'
    vnc.units = 'degree_north'
    vnc[:,:] = grdROMS.lat_psi[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_psi', 'd', ('eta_psi','xi_psi',),zlib=False)
    vnc.long_name = 'Longitude at PSI points'
    vnc.units = 'degree_east'
    vnc[:,:] = grdROMS.lon_psi[startindex:endindex:decimate,startindex:endindex:decimate]
   
    vnc = f1.createVariable('lon_vert', 'd', ('eta_vert','xi_vert',),zlib=False)
    vnc.long_name = 'Longitude at vertices points'
    vnc.units = 'degree_east'
    vnc[:,:] = grdROMS.lon_vert[::decimate,::decimate]
    
    vnc = f1.createVariable('lat_vert', 'd', ('eta_vert','xi_vert',),zlib=False)
    vnc.long_name = 'Latitude at vertices points'
    vnc.units = 'degree_north'
    vnc[:,:] = grdROMS.lat_vert[::decimate,::decimate]

    vnc = f1.createVariable('x_rho', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'X location of RHO points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.x_rho[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('y_rho', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Y location of RHO points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.y_rho[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('x_u', 'd', ('eta_u','xi_u',),zlib=False)
    vnc.long_name = 'X location of U points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.x_u[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('y_u', 'd', ('eta_u','xi_u',),zlib=False)
    vnc.long_name = 'Y location of U points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.y_u[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('x_v', 'd', ('eta_v','xi_v',),zlib=False)
    vnc.long_name = 'X location of V points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.x_v[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('y_v', 'd', ('eta_v','xi_v',),zlib=False)
    vnc.long_name = 'Y location of V points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.y_v[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('x_psi', 'd', ('eta_psi','xi_psi',),zlib=False)
    vnc.long_name = 'X location of PSI points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.x_psi[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('y_psi', 'd', ('eta_psi','xi_psi',),zlib=False)
    vnc.long_name = 'Y location of PSI points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.y_psi[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('x_vert', 'd', ('eta_vert','xi_vert',),zlib=False)
    vnc.long_name = 'X location of VERT points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.x_vert[::decimate,::decimate]

    vnc = f1.createVariable('y_vert', 'd', ('eta_vert','xi_vert',),zlib=False)
    vnc.long_name = 'Y location of VERT points'
    vnc.units = 'meter'
    vnc[:,:] = grdROMS.y_vert[::decimate,::decimate]

    vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Final bathymetry at RHO points'
    vnc.units = 'meter'
    vnc.field = "bath, scalar"
    vnc[:,:] = grdROMS.h[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('hraw', 'd', ('bath','eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Raw bathymetry at RHO-points'
    vnc.units = 'meter'
    vnc.field = "bath, scalar"
    vnc[:,:,:] = grdROMS.hraw[:,startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('f', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'Coriolis parameter at RHO points'
    vnc.units = 'second-1'
    vnc.field = "Coriolis, scalar"
    vnc[:,:] = grdROMS.f[startindex:endindex:decimate,startindex:endindex:decimate]
    
    # Calculate grid metrics
    dndx,dmde,pm,pn = calculateGRDMetrics.calculateGridMetrics(grdROMS,True,decimate,startindex,endindex)
   
    vnc = f1.createVariable('pm', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'curvilinear coordinate metric in XI'
    vnc.units = 'meter-1'
    vnc.field = "pm, scalar"
    vnc[:,:] = pm

    print("Average DX in meters: ",1/np.average(pm,axis=None))
     
    vnc = f1.createVariable('pn', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'curvilinear coordinate metric in ETA'
    vnc.units = 'meter-1'
    vnc.field = "pn, scalar"
    vnc[:,:] = pn

    print("Average DY in meters: ",1/np.average(pn,axis=None))
    
    vnc = f1.createVariable('dmde', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'XI derivative of inverse metric factor pn'
    vnc.units = 'meter'
    vnc.field = "dmde, scalar"
    vnc[:,:] = dmde
    
    vnc = f1.createVariable('dndx', 'd', ('eta_rho','xi_rho',),zlib=False)
    vnc.long_name = 'ETA derivative of inverse metric factor pm'
    vnc.units = 'meter'
    vnc.field = "dndx, scalar"
    vnc[:,:] = dndx
    
    vnc=f1.createVariable('angle','d',('eta_rho','xi_rho'),zlib=False)
    vnc.long_name = "angle between xi axis and east"
    vnc.units = "radian" 
    vnc[:,:]=grdROMS.angle[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho'),zlib=False)
    vnc.long_name = "mask on RHO-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0

    values=grdROMS.mask_rho[startindex:endindex:decimate,startindex:endindex:decimate]
    infile="/Users/trondkr/Projects/KINO/GRID/mask_change.txt"
    fmask=open(infile,"r")

    lines=fmask.readlines()
    import string
    for line in lines:
        l=string.split(line," ")
  
        i=int(float(l[0].strip()))
        j=int(float(l[1].strip()))
        m=float(l[2].strip())
        print("Changing %s %s from %s to %s"%(j,i,values[j,i],m))
        values[j,i]=m
    vnc[:,:]=values[:,:]

    vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u'),zlib=False)
    vnc.long_name = "mask on U-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v'),zlib=False)
    vnc.long_name = "mask on V-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_psi','d',('eta_psi', 'xi_psi'),zlib=False)
    vnc.long_name = "mask on PSI-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_psi[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('s_rho', 'd', ('s_rho',), zlib=False)
    vnc.long_name = "S-coordinate at RHO-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_rho, scalar"
    vnc[:] = grdROMS.s_rho[:]

    vnc = f1.createVariable('s_w', 'd', ('s_w',), zlib=False)
    vnc.long_name = "S-coordinate at W-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_w, scalar"
    vnc[:] = grdROMS.s_w[:]

    vnc = f1.createVariable('Cs_r', 'd', ('s_rho',), zlib=False)
    vnc.long_name = "S-coordinate stretching curves at RHO-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_rho, scalar"
    vnc[:] = grdROMS.Cs_rho[:]

    vnc = f1.createVariable('Cs_w', 'd', ('s_w',), zlib=False)
    vnc.long_name = "S-coordinate stretching curves at W-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_w, scalar"
    vnc[:] = grdROMS.Cs_w[:]

    vnc = f1.createVariable('hc', 'd')
    vnc.long_name = "S-coordinate parameter, critical depth";
    vnc.units = "meter"
    vnc[:] = grdROMS.hc

    vnc = f1.createVariable('xl', 'd')
    vnc.long_name = "domain length in the XI-direction";
    vnc.units = "meter"
    vnc[:] = grdROMS.xl

    vnc = f1.createVariable('el', 'd')
    vnc.long_name = "domain length in the ETA-direction";
    vnc.units = "meter"
    vnc[:] = grdROMS.el

    vnc = f1.createVariable('Tcline', 'd')
    vnc.long_name = "S-coordinate surface/bottom layer width";
    vnc.units = "meter"
    vnc[:] = grdROMS.Tcline

    vnc = f1.createVariable('theta_s', 'd')
    vnc.long_name = "S-coordinate surface control parameter";
    vnc[:] = grdROMS.theta_s
    
    vnc = f1.createVariable('theta_b', 'd')
    vnc.long_name = "S-coordinate bottom control parameter";
    vnc[:] = grdROMS.theta_b

    vnc = f1.createVariable('spherical', 'c')
    vnc.long_name = "Grid type logical switch";
    vnc.flag_values = "T, F" ;
    vnc.flag_meanings = "spherical Cartesian" ;

    vnc[:] = grdROMS.spherical


    f1.close()

    print("Creating new decimated grid file: %s"%(outfilename))  
    



