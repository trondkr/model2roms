"""
This class creates an object based on the input file structure. Currently the class takes
two types of structures as input: SODA or ROMS
"""
import logging
from datetime import datetime

import numpy as np
import xarray as xr

import IOverticalGrid

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2008, 12, 9)
__modified__ = datetime(2021, 3, 23)
__version__ = "1.5"
__status__ = "Development"


class Grd:

    def __init__(self, grdtype, confM2R):
        """
        The object is initialised and created through the __init__ method
        As an example of how to use, these lines return a grid object called grdTEST:
        => import grd
        => grdTEST = grd.grdClass("grdfilename","ROMS")
        """
        self.type = grdtype
        self.grdName = confM2R.outgrid_name
        self.realm = confM2R.realm

        logging.info("[M2R_grd] Creating init for grid object {}".format(confM2R.outgrid_name))
        logging.info("[M2R_grd]---> Initialized GRD object for grid type {}\n".format(self.type))

    def create_object(self, confM2R, grd_filename):
        """
        This method creates a new object by reading the grd input file. All
        dimensions (eta, xi, lon, lat etc.) are defined here and used throughout these scripts.
        Also, the depth matrix is calculated in this function by calling IOverticalGrid.py (ROMS grid only). For
        input model depths, the depth array is a one dimensional. If input data has 2 or 3 dimensions, this
        has to be accounted for throughout the model2roms package as one dimension is currently only supported.

        Object types currently supported: STATION, GLORYS, SODA3

        Trond Kristiansen, 18.11.2009, edited 03.01.2017, 23.03.2021
        """
        ds = xr.open_dataset(grd_filename)

        if self.type == 'FORCINGDATA':

            logging.info("[M2R_grd] ---> Assuming {} grid type for {}".format(confM2R.grd_type, self.type))
            logging.info("[M2R_grd] ---> Using dimension names {} and {} and {}".format(confM2R.lon_name,
                                                                                        confM2R.lat_name,
                                                                                        confM2R.depth_name))

            self.lon = ds[str(confM2R.lon_name)][:]
            self.lat = ds[str(confM2R.lat_name)][:]
            self.h = ds[str(confM2R.depth_name)][:]
            self.nlevels = len(self.h)
            self.fillval = -9.99e+33
            self.hc = None

            if self.lon.ndim == 1:
                self.lon, self.lat = np.meshgrid(self.lon, self.lat)

            # Create grid for ESMF interpolation

            self.esmfgrid = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                      is_sphere=True, coord_names=[str(confM2R.lon_name), str(confM2R.lat_name)],
                                      add_mask=False)
            self.esmfgrid_u = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                        is_sphere=True,
                                        coord_names=[str(confM2R.lon_name_u), str(confM2R.lat_name_u)],
                                        add_mask=False)
            self.esmfgrid_v = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                        is_sphere=True,
                                        coord_names=[str(confM2R.lon_name_v), str(confM2R.lat_name_v)],
                                        add_mask=False)

            if confM2R.ocean_indata_type == 'SODA3':
                self.fillval = -1.e+20
            if confM2R.ocean_indata_type == 'SODA3_5DAY':
                self.fillval = -1.e+20
            if confM2R.ocean_indata_type == 'GLORYS':
                self.fillval = 9.96921e+36

            if confM2R.ocean_indata_type == 'NORESM':
                # self.h = ds["depth"][:]
                self.h = np.asarray([0, 5, 10, 15, 20, 25, 30, 40, 50, 62.5, 75, 87.5, 100, 112.5, 125,
                                     137.5, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600,
                                     650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
                                     1300, 1350, 1400, 1450, 1500, 1625, 1750, 1875, 2000, 2250, 2500, 2750,
                                     3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750,
                                     6000, 6250, 6500, 6750])
                self.fillval = 32768
                self.nlevels = len(self.h)

            IOverticalGrid.get_z_levels(self)

        if self.type == 'STATION':
            self.lon = ds[confM2R.lon_name][:]
            self.lat = ds[confM2R.lat_name][:]
            self.h = ds[confM2R.depth_name][:]
            self.time = ds[confM2R.time_name][:]

            self.Lp = 1
            self.Mp = 1
            self.fillval = -9.99e+33

        if self.type in ['ROMS']:

            self.write_clim = True
            self.write_bry = True
            self.write_init = True
            self.write_stations = False

            self.lonname = 'lon_rho'
            self.latname = 'lat_rho'

            """
            Set initTime to 1 if you dont want the first time-step to be
            the initial field (no ubar and vbar if time=0)
            """

            self.inittime = 0
            self.ocean_time = 0
            self.NT = 2
            self.tracer = self.NT

            self.message = None  # Used to store the date for printing to screen (IOwrite.py)
            self.time = 0
            self.reftime = 0
            self.grdtype = 'regular'
            self.mask_rho = ds["mask_rho"][:, :]
            self.lon_rho = ds["lon_rho"][:, :]
            self.lat_rho = ds["lat_rho"][:, :]
            self.h = ds["h"][:, :]

            masked_h = np.where(self.h > 0, self.h, self.h.max())

            self.hmin = masked_h.min()
            if "Vtransform" in ds.variables:
                self.vtransform = ds["Vtransform"].values
            else:
                self.vtransform = confM2R.vtransform

            if "s_rho" in ds.variables:
                self.s_rho = ds["s_rho"].values
                self.nlevels = len(self.s_rho)
            else:
                self.nlevels = confM2R.nlevels

            if "Vstretching" in ds.variables:
                self.vstretching = ds["Vstretching"].values
            if "theta_s" in ds.variables:
                self.theta_s = ds["theta_s"].values
            else:
                self.theta_s = confM2R.theta_s
            if "theta_b" in ds.variables:
                self.theta_b = ds["theta_b"].values
            else:
                self.theta_b = confM2R.theta_b
            if "Tcline" in ds.variables:
                self.tcline = ds["Tcline"].values
            else:
                self.tcline = confM2R.tcline
            if "hc" in ds.variables:
                self.hc = ds["hc"].values
            else:
                self.hc = confM2R.hc

            if self.vtransform == 1:
                self.hc = min(self.hmin, self.tcline)
                self.hc = self.tcline
                if self.tcline > self.hmin:
                    print('Vertical transformation parameters are not defined correctly in either gridid.txt '
                          'or in the history files: \n Tc\
                           line = %d and hmin = %d. \n You need to make sure that '
                          'tcline <= hmin when using transformation 1.' % (
                              self.tcline, self.hmin))
            else:
                self.hc = self.tcline

            zeta = None
            if zeta is None:
                self.zeta = np.zeros(self.h.shape)
            else:
                self.zeta = zeta

            # for findvar in ds:
            #    if findvar=="hraw":
            #        self.hraw     = ds["hraw"][:,:,:]

            self.lon_u = ds["lon_u"][:, :]
            self.lat_u = ds["lat_u"][:, :]
            self.mask_u = ds["mask_u"][:, :]
            for findvar in ds:
                if findvar == "lon_vert":
                    self.lon_vert = ds["lon_vert"][:, :]
                    self.lat_vert = ds["lat_vert"][:, :]

            for findvar in ds:
                if findvar == "x_rho":
                    self.x_rho = ds["x_rho"][:, :]
                    self.y_rho = ds["y_rho"][:, :]

            for findvar in ds:
                if findvar == "x_u":
                    self.x_u = ds["x_u"][:, :]
                    self.y_u = ds["y_u"][:, :]

            for findvar in ds:
                if findvar == "x_v":
                    self.x_v = ds["x_v"][:, :]
                    self.y_v = ds["y_v"][:, :]

            for findvar in ds:
                if findvar == "x_psi":
                    self.x_psi = ds["x_psi"][:, :]
                    self.y_psi = ds["y_psi"][:, :]

            for findvar in ds:
                if findvar == "x_vert":
                    self.x_vert = ds["x_vert"][:, :]
                    self.y_vert = ds["y_vert"][:, :]

            for findvar in ds:
                if findvar == "xl":
                    self.xl = ds["xl"]
                    self.el = ds["el"]

            for findvar in ds:
                if findvar == "dmde":
                    self.dmde = ds["dmde"][:, :]
                    self.dndx = ds["dndx"][:, :]

            self.lon_v = ds["lon_v"][:, :]
            self.lat_v = ds["lat_v"][:, :]
            self.mask_v = ds["mask_v"][:, :]

            #   self.spherical = ds["spherical"][:]

            self.lon_psi = self.lon_u[:-1, :]
            self.lat_psi = self.lat_v[:, :-1]
            self.mask_psi = self.mask_v[:, :-1]

            # self.f = ds["f"][:, :]
            self.angle = ds["angle"][:, :]

            self.pm = ds["pm"][:, :]
            self.invpm = 1.0 / np.asarray(ds["pm"][:, :])
            self.pn = ds["pn"][:, :]
            self.invpn = 1.0 / np.asarray(ds["pn"][:, :])

            self.Lp = len(self.lat_rho[1, :])
            self.Mp = len(self.lat_rho[:, 1])

            self.fillval = -9.99e33

            self.eta_rho = self.Mp
            self.eta_u = self.Mp
            self.eta_v = self.Mp - 1
            self.eta_psi = self.Mp - 1
            self.xi_rho = self.Lp
            self.xi_u = self.Lp - 1
            self.xi_v = self.Lp
            self.xi_psi = self.Lp - 1

            # Boolean to check if we need to initialize the CLIM file before writing
            self.ioClimInitialized = False
            self.ioInitInitialized = False

            if self.lon_rho.ndim == 1:
                self.lon_rho, self.lat_rho = np.meshgrid(self.lon_rho, self.lat_rho)
                self.lon_u, self.lat_u = np.meshgrid(self.lon_u, self.lat_u)
                self.lon_v, self.lat_v = np.meshgrid(self.lon_v, self.lat_v)

            # Setup the vertical coordinate system
            IOverticalGrid.calculateVgrid(self)

            self.esmfgrid_u = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                        coord_names=['lon_u', 'lat_u'],
                                        is_sphere=True,
                                        add_mask=False)
            self.esmfgrid_v = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                        is_sphere=True,
                                        coord_names=['lon_v', 'lat_v'],
                                        add_mask=False)
            self.esmfgrid = ESMF.Grid(filename=grd_filename, filetype=ESMF.FileFormat.GRIDSPEC,
                                      is_sphere=True,
                                      coord_names=[self.lonname, self.latname],
                                      add_mask=False)

    def getdims(self):
        if self.type in ["ROMS"]:
            self.Lp = len(self.lat_rho[1, :])
            self.Mp = len(self.lat_rho[:, 1])

        if self.type in ['FORCINGDATA']:
            self.Lp = len(self.lat[1, :])
            self.Mp = len(self.lat[:, 1])

        self.M = self.Mp - 1
        self.L = self.Lp - 1
