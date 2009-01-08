

"""
This class creates an object based on the input file structure. Currently the class takes
two types of structures as input: SODA or ROMS
"""

import os, sys
from datetime import datetime
from netCDF4 import Dataset
import numpy as np

import IOverticalGrid
import printObject

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 12, 9)
__modified__ = datetime(2008, 12, 9)
__version__  = "1.1"
__status__   = "Development"


class grdClass:

    def __init__(self,grdfilename,type,log):
        """
        The object is initialised and created through the __init__ method
        """
        self.grdfilename= grdfilename
      
        self.log=log
        self.type=type
        
        self.openNetCDF()
        self.createObject()
        self._getDims()
        
        print 'Generating GRD object for grid type %s'%(self.type)
        
        if self.type=="ROMS":
            self.Nlevels=30
            self.theta_s=5.0
            self.theta_b=0.9
            self.Tcline=50.0
            self.hc=10.0
            self.ocean_time=1
            # Number of tracers
            self.NT=2
            self.tracer=self.NT
    
  
    def openNetCDF(self):
        """
        Open the netCDF file and store the contents in arrays associated with variable names
        """
        try:
            self.cdf = Dataset(self.grdfilename,"r")
        except IOError:
            print 'Could not open file %s'%(self.grdfilename)
            print 'Exception caught in: openNetCDF(grdfilename, log)'
        if self.log is True:  
            print '\n---> Opened input file %s'%(self.grdfilename)
        
  
    def createObject(self):
        """
        This method creates a new object by reading the grd input file
        """
        if self.type=='SODA':
            self.lon = self.cdf.variables["LON"][:]
            self.lat = self.cdf.variables["LAT"][:]
            self.depth = self.cdf.variables["DEPTH"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33
            
            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)
           
            
        if self.type=='ROMS':
            self.lon_rho  = self.cdf.variables["lon_rho"][:]
            self.lat_rho  = self.cdf.variables["lat_rho"][:]
            self.depth    = self.cdf.variables["h"][:,:]
            self.mask_rho = self.cdf.variables["mask_rho"][:,:]
          
            self.lon_u  = self.cdf.variables["lon_u"][:]
            self.lat_u  = self.cdf.variables["lat_u"][:]
            self.mask_u = self.cdf.variables["mask_u"][:,:]
            
            self.lon_v  = self.cdf.variables["lon_v"][:]
            self.lat_v  = self.cdf.variables["lat_v"][:]
            self.mask_v = self.cdf.variables["mask_v"][:,:]
            
            self.lon_psi  = self.cdf.variables["lon_psi"][:]
            self.lat_psi  = self.cdf.variables["lat_psi"][:]
            self.mask_psi = self.cdf.variables["mask_psi"][:,:]
            
            self.f  = self.cdf.variables["f"][:]
            self.xl  = self.cdf.variables["xl"][:]
            self.el  = self.cdf.variables["el"][:]
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
            
            if np.rank(self.lon_rho)==1:
                    self.lon_rho, self.lat_rho = np.meshgrid(self.lon_rho,self.lat_rho)
                    self.lon_u, self.lat_u = np.meshgrid(self.lon_u,self.lat_u)
                    self.lon_v, self.lat_v = np.meshgrid(self.lon_v,self.lat_v)
    
 
    def _getDims(self):
        if self.type=="ROMS":
            self.Lp=len(self.lat_rho[1,:])
            self.Mp=len(self.lat_rho[:,1])
        if self.type=="SODA":
            self.Lp=len(self.lat[1,:])
            self.Mp=len(self.lat[:,1])
        self.M =self.Mp-1
        self.M =self.Mp-1
        self.L =self.Lp-1
        
    #def rotateUV(self):
    #    angle = self.angle;
    #    
    #    self.u =  u*cos(angle) + v*sin(angle);
    #    vout = -u.*sin(angle) + v*cos(angle);
