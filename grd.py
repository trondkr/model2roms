

"""
This class creates an object based on the input file structure. Currently the class takes
two types of structures as input: SODA or ROMS
"""

from datetime import datetime
from netCDF4 import Dataset
import numpy as np

import IOverticalGrid
try:
    import ESMF
except ImportError:
    print "Could not find module ESMF"
    pass

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 12, 9)
__modified__ = datetime(2014, 10, 23)
__version__  = "1.2"
__status__   = "Development"


class grdClass:

    def __init__(self,grdfilename,type,name,useESMF):
        """
        The object is initialised and created through the __init__ method
        As an example of how to use, these lines return a grid object called grdTEST:
        => import grd
        => grdTEST = grd.grdClass("grdfilename","ROMS")
        """

        self.grdfilename= grdfilename
        self.type=type
        self.grdName=name
        self.useESMF=useESMF

        self.openNetCDF()
        self.createObject()
        self.getDims()

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
        if self.type=='AVERAGE':
            # Average file is a self made climatology file foe the North Atlantic (created with averageSoda.py)
            # and used in individual-based modeling (ibm.py)
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["lon"][:]
            self.lat = self.cdf.variables["lat"][:]
            self.lonName='lon'
            self.latName='lat'
            self.depth = self.cdf.variables["depth"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33

            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)

        if self.type=='WOAMONTHLY':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["lon"][:]
            self.lat = self.cdf.variables["lat"][:]
            self.lonName='lon'
            self.latName='lat'
            self.depth = self.cdf.variables["depth"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=9.96921e+36

            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)

            IOverticalGrid.get_z_levels(self)

        if self.type=='SODA':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["LON"][:]
            self.lat = self.cdf.variables["LAT"][:]
            self.lonName='LON'
            self.latName='LAT'
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
            self.lonName='lon'
            self.latName='lat'
            self.depth = self.cdf.variables["depth"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=-9.99e+33

            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)

            IOverticalGrid.get_z_levels(self)
            
        if self.type=='GLORYS':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["nav_lon"][:]
            self.lat = self.cdf.variables["nav_lat"][:]
            self.lonName='nav_lon'
            self.latName='nav_lat'
            #NOTE spelling error for depth in netcdf files
            self.depth = self.cdf.variables["deptht"][:]
            self.Nlevels = len(self.depth)
            self.fill_value=9.96921e+36

            if np.rank(self.lon)==1:
                    self.lon, self.lat = np.meshgrid(self.lon,self.lat)

            IOverticalGrid.get_z_levels(self)

        if self.type=='NORESM':
            self.grdType  = 'regular'
            print '---> Assuming %s grid type for %s'%(self.grdType,self.type)
            self.lon = self.cdf.variables["plon"][:]
            self.lat = self.cdf.variables["plat"][:]
            self.lonName='plon'
            self.latName='plat'
            self.fieldSrc='NaN'
            self.fieldDst_rho='NaN'
            self.fieldDst_u='NaN'
            self.fieldDst_v='NaN'

           # self.depth = self.cdf.variables["depth"][:]
            self.depth=np.asarray([0, 5, 10, 15, 20, 25, 30, 40, 50, 62.5, 75, 87.5, 100, 112.5, 125,
    137.5, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600,
    650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
    1300, 1350, 1400, 1450, 1500, 1625, 1750, 1875, 2000, 2250, 2500, 2750,
    3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750,
    6000, 6250, 6500, 6750])


            self.Nlevels = len(self.depth)
            self.fill_value=-32768

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

        if self.type=='MARMAP':

            self.lon = 0.0
            self.lat = 0.0
            self.depth = 0.0
            self.time = 0.0

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
            self.theta_s=5.0
            self.theta_b=0.4
            self.Tcline=20.0
            self.hc=20.0
            self.vars=[]
            self.lonName='lon_rho'
            self.latName='lat_rho'

            """Set initTime to 1 if you dont want the first timestep to be
            the initial field (no ubar and vbar if time=0)"""

            self.initTime=0
            self.ocean_time=0
            self.NT=2
            self.tracer=self.NT

            """Parameters that are used to fill in gaps in interpolated fields
            due to mismatch of iput and output mask (used in cleanArray.f90)"""

            self.maxval=1000
            self.minDistPoints=5
            self.maxDistHorisontal=100
            self.maxDistVertical=100

            self.message  = None  # Used to store the date for printing to screen (IOwrite.py)
            self.time     = 0
            self.reftime  = 0
            self.grdType  = 'regular'
            self.lon_rho  = self.cdf.variables["lon_rho"][:,:]
            self.lat_rho  = self.cdf.variables["lat_rho"][:,:]

            self.depth    = self.cdf.variables["h"][:,:]

            self.mask_rho = self.cdf.variables["mask_rho"][:,:]

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

            if (self.useESMF):

                self.esmfgrid_u = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                      is_sphere=True, coord_names=['lon_u','lat_u'], add_mask=False)
                self.esmfgrid_v = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                      is_sphere=True, coord_names=['lon_v','lat_v'], add_mask=False)

        # Create grid for ESMF interpolation
        if (self.useESMF):
            self.esmfgrid = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                      is_sphere=True, coord_names=[self.lonName, self.latName], add_mask=False)



    def getDims(self):
        if self.type=="ROMS":
            self.Lp=len(self.lat_rho[1,:])
            self.Mp=len(self.lat_rho[:,1])
        if self.type in ['SODA','SODAMONTHLY','AVERAGE','GLORYS','NORESM','WOAMONTHLY']:
            self.Lp=len(self.lat[1,:])
            self.Mp=len(self.lat[:,1])

        self.M =self.Mp-1
        self.L =self.Lp-1
