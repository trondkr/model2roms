

"""
This class creates an object based on the input file structure. Currently the class takes
two types of structures as input: SODA or ROMS
"""

import os, sys
from datetime import datetime
from netCDF4 import Dataset
import numpy as np

import IOverticalGrid

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 12, 9)
__modified__ = datetime(2009, 3, 10)
__version__  = "1.1"
__status__   = "Development"


class grdClass:

    def __init__(self,grdfilename,type):
        """
        The object is initialised and created through the __init__ method
        As an example of how to use, these lines return a grid object called grdTEST:
        => import grd
        => grdTEST = grd.grdClass("grdfilename","ROMS")
        """
        
        self.grdfilename= grdfilename
        self.type=type
        
        self.openNetCDF()
        self.createObject()
        self.getDims()
        
        if grdfilename=="/Users/trond/Projects/arcwarm/nordic/AA_10km_grid.nc":
            self.grdName='NA'
        elif grdfilename=="/Users/trond/Projects/Roms/Nordic/Inputs/imr_nordic_4km.nc":
            self.grdName='Nordic'
        elif grdfilename=='/Users/trond/Projects/Nathan/NoMed47_GRID_Global.nc':
            self.grdName='NA_Nathan'
        elif grdfilename=='/Users/trond/Projects/Nathan/GOM_GRID_Global.nc':
            self.grdName='GOM_Nathan'
        elif grdfilename=='/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc':
            self.grdName='Nordic2'
        elif grdfilename=='/Users/trond/Projects/Roms/GOMsmall/Inputs/ns8_grd.nc':
            self.grdName='NorthSea'
        elif grdfilename=='/Users/trond/Projects/Roms/Julia/NATLC/Data/natl_40km.nc':
            self.grdName='NATL'
            
        else:
            self.grdName='GOM'
            
        print '---> Generated GRD object for grid type %s'%(self.type)
    
  
    def openNetCDF(self):
        """Open the netCDF file and store the contents in arrays associated with variable names"""
        try:
            self.cdf = Dataset(self.grdfilename,"r")
        except IOError:
            print 'Could not open file %s'%(self.grdfilename)
            print 'Exception caught in: openNetCDF(grdfilename)'
        
  
    def createObject(self):
        """
        This method creates a new object by reading the grd input file. All
        dimensions (eta, xi, lon, lat etc.) are defined here and used througout these scripts.
        Also, the depth matrix is calculated in this function by calling IOverticalGrid.py (ROMS grid only). For
        input model depths, the depth array is a one dimensional. If input data has 2 or 3 dimensions, this
        has to be accounted for througout the soda2roms package as one dimension is currently only supported.
        
        Object types currently supported: HYCOM,SODA,SODAMONTHLY,STATION,ROMS
        
        Trond Kristiansen, 18.11.2009.
        """
        if self.type=='SODA':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["LON"][:]
            self.lat = self.cdf.variables["LAT"][:]
            self.depth = self.cdf.variables["DEPTH"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33
            
            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)
            
            IOverticalGrid.get_z_levels(self)
            
        if self.type=='SODAMONTHLY':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["lon"][:]
            self.lat = self.cdf.variables["lat"][:]
            self.depth = self.cdf.variables["depth"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33
            
            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)
            
            IOverticalGrid.get_z_levels(self)
        
        if self.type=='HYCOM':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["Longitude"][:]
           
            if self.lon.max() > 360:
                if np.rank(self.lon)==1:
                    self.lon[:]=np.where(self.lon>360,self.lon[:]-360,self.lon[:])
                else:
                    self.lon[:,:]=np.where(self.lon>360,self.lon[:,:]-360,self.lon[:,:])
            
            
            self.lat = self.cdf.variables["Latitude"][:]
            self.depth = self.cdf.variables["Depth"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33
            
            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)
            
            IOverticalGrid.get_z_levels(self)
            
        if self.type=='STATION':
            
            self.lon = self.cdf.variables["lon"][:]
            self.lat = self.cdf.variables["lat"][:]
            self.depth = self.cdf.variables["depth"][:]
            self.time = self.cdf.variables["time"][:]
            
            self.Lp=1
            self.Mp=1
            self.fill_value=-9.99e+33
            
        if self.type=='ROMS':
            
            self.write_clim=True
            self.write_bry=True
            self.write_init=True
            self.write_stations=False
    
            """Define varibales needed to set the apropriate vertical stretching parameters. Note
            that according to the ROMS forum
            (https://www.myroms.org/forum/viewtopic.php?f=23&t=1254&hilit=critical+depth+tcline&sid=ec98a9e63e7857e2615b9182af752cde)
            the value of Tcline should now be equal to hc"""
            
            self.vstretching=2
            self.Nlevels=35
            self.theta_s=1.0
            self.theta_b=3.0
            self.Tcline=10.0
            self.hc=10.0
            self.vars=[]
            
            """Set initTime to 1 if you dont want the first timestep to be
            the initial field (no ubar and vbar if time=0)"""
            
            self.initTime=1 
            self.ocean_time=1
            self.NT=2
            self.tracer=self.NT
            
            """Parameters that are used to fill in gaps in interpolated fields
            due to mismatch of iput and output mask (used in cleanArray.f90)"""
            
            self.maxval=1000
            self.minDistPoints=5
            self.maxDistHorisontal=50
            self.maxDistVertical=50
            
            self.message  = None  # Used to store the date for printing to screen (IOwrite.py)
            self.time     = 0
            self.reftime  = 0
            self.grdType  = 'regular'
            self.lon_rho  = self.cdf.variables["lon_rho"][:,:]
            self.lat_rho  = self.cdf.variables["lat_rho"][:,:]
            self.depth    = self.cdf.variables["h"][:,:]
          
            self.mask_rho = self.cdf.variables["mask_rho"][:,:]
            # Nathan fixes
            # Need to convert all longitude values in grid that are less than 360 and larger than 180 to negative values.
            #self.lon_rho=self.lon_rho-360
            # Cannot have undefined values in depth matrix. This messes up the z_r calculautions. Therefore,
            # convert all depth values that are not valid (e.g. >10000 m) to 0 and set the shallowest depth to self.hc.
            # Do this in the gird file and remove from here.
            
            #self.depth[:,:]=np.where(self.depth[:,:]>9000,self.hc,self.depth[:,:])    
            #self.depth[:,:]=np.where(self.depth[:,:]<self.hc,self.hc,self.depth[:,:])
            
           
            self.lon_u  = self.cdf.variables["lon_u"][:,:]
            self.lat_u  = self.cdf.variables["lat_u"][:,:]
            self.mask_u = self.cdf.variables["mask_u"][:,:]
            
            self.lon_v  = self.cdf.variables["lon_v"][:,:]
            self.lat_v  = self.cdf.variables["lat_v"][:,:]
            self.mask_v = self.cdf.variables["mask_v"][:,:]
            
            self.lon_psi  = self.lon_u[:-1,:]
            self.lat_psi  = self.lat_v[:,:-1]
            self.mask_psi = self.mask_v[:,:-1]
            
            self.f  = self.cdf.variables["f"][:,:]
            self.angle  = self.cdf.variables["angle"][:,:]
            
            self.pm  = self.cdf.variables["pm"][:,:]
            self.invpm  = 1.0/np.asarray(self.cdf.variables["pm"][:,:])
            self.pn  = self.cdf.variables["pn"][:,:]
            self.invpn  = 1.0/np.asarray(self.cdf.variables["pn"][:,:])
            
            self.Lp=len(self.lat_rho[1,:])
            self.Mp=len(self.lat_rho[:,1])
            
            self.fill_value=-9.99e+33
            
            self.eta_rho = self.Mp
            self.eta_u   = self.Mp
            self.eta_v   = self.Mp-1
            self.eta_psi   = self.Mp-1
            self.xi_rho  = self.Lp
            self.xi_u    = self.Lp-1
            self.xi_v    = self.Lp
            self.xi_psi    = self.Lp-1
            
            """Boolean to check if we need to initialize the CLIM file before writing"""
            self.ioClimInitialized=False
            self.ioInitInitialized=False
            
            if np.rank(self.lon_rho)==1:
                    self.lon_rho, self.lat_rho = np.meshgrid(self.lon_rho,self.lat_rho)
                    self.lon_u, self.lat_u = np.meshgrid(self.lon_u,self.lat_u)
                    self.lon_v, self.lat_v = np.meshgrid(self.lon_v,self.lat_v)
                    
            IOverticalGrid.calculate_z_r(self)
            IOverticalGrid.calculate_z_w(self)
           
        
    def getDims(self):
        if self.type=="ROMS":
            self.Lp=len(self.lat_rho[1,:])
            self.Mp=len(self.lat_rho[:,1])
        if self.type=="SODA" or self.type=='SODAMONTHLY':
            self.Lp=len(self.lat[1,:])
            self.Mp=len(self.lat[:,1])
        if self.type=="HYCOM":
            self.Lp=len(self.lat[1,:])
            self.Mp=len(self.lat[:,1])
        self.M =self.Mp-1
        self.L =self.Lp-1
        
