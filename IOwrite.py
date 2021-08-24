import os
import time
from datetime import datetime

import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2009, 3, 2)
__modified__ = datetime(2014, 10, 23)
__version__ = "1.1"
__status__ = "Development"


def help():
    """
    This function generates a CLIM file from scratch. The variables created include:
    salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is ocean_time which is days
    since 1948/1/1.

    This file is netcdf CF compliant and to check the CLIM file for CF compliancy:
    http://titania.badc.rl.ac.uk/cgi-bin/cf-checker.pl?cfversion=1.0

    (also see: https://www.myroms.org/forum/viewtopic.php?f=30&t=1450&p=5209&hilit=cf+compliant#p5209)

    This function is called from model2roms.py.
    """


def write_clim_file(confM2R, ntime, myvar, data1=None, data2=None, data3=None, data4=None):
    if confM2R.output_format == 'NETCDF4':
        myzlib = True
    else:
        myzlib = False

    grdROMS = confM2R.grdROMS

    if confM2R.grdROMS.ioClimInitialized is False:
        confM2R.grdROMS.ioClimInitialized = True
        if os.path.exists(confM2R.clim_name):
            os.remove(confM2R.clim_name)

        f1 = Dataset(confM2R.clim_name, mode='w', format=confM2R.output_format)
        f1.title = "Climatology forcing file (CLIM) used for forcing the ROMS model"
        f1.description = "Created for grid file: %s" % (confM2R.roms_grid_path)
        f1.grd_file = "Gridfile: %s" % (confM2R.roms_grid_path)
        f1.history = "Created " + time.ctime(time.time())
        f1.source = "{} ({})".format(confM2R.author_name, confM2R.author_email)
        f1.type = "File in %s format created using MODEL2ROMS" % (confM2R.output_format)
        f1.link = "https://github.com/trondkr/model2roms"
        f1.Conventions = "CF-1.0"

        # Define dimensions
        f1.createDimension('xi_rho', grdROMS.xi_rho)
        f1.createDimension('eta_rho', grdROMS.eta_rho)
        f1.createDimension('xi_u', grdROMS.xi_u)
        f1.createDimension('eta_u', grdROMS.eta_u)
        f1.createDimension('xi_v', grdROMS.xi_v)
        f1.createDimension('eta_v', grdROMS.eta_v)
        f1.createDimension('xi_psi', grdROMS.xi_psi)
        f1.createDimension('eta_psi', grdROMS.eta_psi)
        f1.createDimension('s_rho', len(grdROMS.s_rho))
        f1.createDimension('s_w', len(grdROMS.s_w))
        if confM2R.isclimatology:
            f1.createDimension('clim_time', 12)
        else:
            f1.createDimension('ocean_time', None)

        vnc = f1.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Longitude of RHO-points'
        vnc.units = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :] = grdROMS.lon_rho

        vnc = f1.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Latitude of RHO-points'
        vnc.units = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :] = grdROMS.lat_rho

        vnc = f1.createVariable('lon_u', 'd', ('eta_u', 'xi_u',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Longitude of U-points'
        vnc.units = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :] = grdROMS.lon_u

        vnc = f1.createVariable('lat_u', 'd', ('eta_u', 'xi_u',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Latitude of U-points'
        vnc.units = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :] = grdROMS.lat_u

        vnc = f1.createVariable('lon_v', 'd', ('eta_v', 'xi_v',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Longitude of V-points'
        vnc.units = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :] = grdROMS.lon_v

        vnc = f1.createVariable('lat_v', 'd', ('eta_v', 'xi_v',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Latitude of V-points'
        vnc.units = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :] = grdROMS.lat_v

        vnc = f1.createVariable('lat_psi', 'd', ('eta_psi', 'xi_psi',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Latitude of PSI-points'
        vnc.units = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :] = grdROMS.lat_psi

        vnc = f1.createVariable('lon_psi', 'd', ('eta_psi', 'xi_psi',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Longitude of PSI-points'
        vnc.units = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :] = grdROMS.lon_psi

        vnc = f1.createVariable('h', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Bathymetry at RHO-points'
        vnc.units = 'meter'
        vnc.field = "bath, scalar"
        vnc[:, :] = grdROMS.h

        #     vnc = f1.createVariable('f', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        #     vnc.long_name = 'Coriolis parameter at RHO-points'
        #     vnc.units = 'second-1'
        #     vnc.field = "Coriolis, scalar"
        #     vnc[:, :] = grdROMS.f

        vnc = f1.createVariable('pm', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'curvilinear coordinate metric in XI'
        vnc.units = 'meter-1'
        vnc.field = "pm, scalar"
        vnc[:, :] = grdROMS.pm

        vnc = f1.createVariable('pn', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'curvilinear coordinate metric in ETA'
        vnc.units = 'meter-1'
        vnc.field = "pn, scalar"
        vnc[:, :] = grdROMS.pn

        vnc = f1.createVariable('s_rho', 'd', ('s_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate at RHO-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        if grdROMS.vtransform == 2:
            vnc.standard_name = "ocean_s_coordinate_g2"
            vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
        if grdROMS.vtransform == 1:
            vnc.standard_name = "ocean_s_coordinate_g1"
            vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
        vnc.field = "s_rho, scalar"
        vnc[:] = grdROMS.s_rho

        vnc = f1.createVariable('s_w', 'd', ('s_w',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate at W-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        if grdROMS.vtransform == 2:
            vnc.standard_name = "ocean_s_coordinate_g2"
            vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
        if grdROMS.vtransform == 1:
            vnc.standard_name = "ocean_s_coordinate_g1"
            vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
        vnc.field = "s_w, scalar"
        vnc[:] = grdROMS.s_w

        vnc = f1.createVariable('Cs_r', 'd', ('s_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate stretching curves at RHO-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        vnc.field = "s_rho, scalar"
        vnc[:] = grdROMS.Cs_rho

        vnc = f1.createVariable('Cs_w', 'd', ('s_w',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate stretching curves at W-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        vnc.field = "s_w, scalar"
        vnc[:] = grdROMS.Cs_w

        vnc = f1.createVariable('hc', 'd')
        vnc.long_name = "S-coordinate parameter, critical depth";
        vnc.units = "meter"
        vnc[:] = grdROMS.hc

        vnc = f1.createVariable('z_r', 'd', ('s_rho', 'eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "Sigma layer to depth matrix";
        vnc.units = "meter"
        vnc[:, :, :] = grdROMS.z_r

        vnc = f1.createVariable('z_w', 'd', ('s_w', 'eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "Sigma layer to depth matrix";
        vnc.units = "meter"
        vnc[:, :, :] = grdROMS.z_w

        vnc = f1.createVariable('Tcline', 'd')
        vnc.long_name = "S-coordinate surface/bottom layer width"
        vnc.units = "meter"
        vnc[:] = grdROMS.tcline

        vnc = f1.createVariable('theta_s', 'd')
        vnc.long_name = "S-coordinate surface control parameter"
        vnc[:] = grdROMS.theta_s

        vnc = f1.createVariable('theta_b', 'd')
        vnc.long_name = "S-coordinate bottom control parameter"
        vnc[:] = grdROMS.theta_b

        vnc = f1.createVariable('angle', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "angle between xi axis and east"
        vnc.units = "radian"
        vnc[:, :] = grdROMS.angle

        # Now start creating variables for regular climatology/bry/init creations
        if not confM2R.isclimatology:
            v_time = f1.createVariable('ocean_time', 'd', ('ocean_time',), zlib=myzlib, fill_value=grdROMS.fillval)
            v_time.long_name = 'seconds since 1948-01-01 00:00:00'
            v_time.units = 'seconds since 1948-01-01 00:00:00'
            v_time.field = 'time, scalar, series'
            if (confM2R.ocean_indata_type == "NORESM"):
                v_time.calendar = 'noleap'
            else:
                v_time.calendar = 'standard'

            v_u = f1.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
            v_u.long_name = "u-momentum component"
            v_u.units = "meter second-1"
            v_u.time = "ocean_time"
            v_u.field = "u-velocity, scalar, series"
            #v_u.missing_value = grdROMS.fillval

            v_v = f1.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
            v_v.long_name = "v-momentum component"
            v_v.units = "meter second-1"
            v_v.time = "ocean_time"
            v_v.field = "v-velocity, scalar, series"
            #v_v.missing_value = grdROMS.fillval

            v_salt = f1.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_salt.long_name = "salinity"
            v_salt.time = "ocean_time"
            v_salt.field = "salinity, scalar, series"
            #v_salt.missing_value = grdROMS.fillval

            v_temp = f1.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_temp.long_name = "potential temperature"
            v_temp.units = "Celsius"
            v_temp.time = "ocean_time"
            v_temp.field = "temperature, scalar, series"
            #v_temp.missing_value = grdROMS.fillval

            v_ssh = f1.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
            v_ssh.long_name = "sea level"
            v_ssh.units = "meter"
            v_ssh.time = "ocean_time"
            v_ssh.field = "sea level, scalar, series"
            #v_ssh.missing_value = grdROMS.fillval

            v_ubar = f1.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_ubar.long_name = "u-2D momentum"
            v_ubar.units = "meter second-1"
            v_ubar.time = "ocean_time"
            v_ubar.field = "u2-D velocity, scalar, series"
            #v_ubar.missing_value = grdROMS.fillval

            v_vbar = f1.createVariable('vbar', 'f', ('ocean_time', 'eta_v', 'xi_v',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_vbar.long_name = "v-2D momentum"
            v_vbar.units = "meter second-1"
            v_vbar.time = "ocean_time"
            v_vbar.field = "v2-D velocity, scalar, series"
            #v_vbar.missing_value = grdROMS.fillval

            if confM2R.write_ice:
                ageice = f1.createVariable('ageice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                           fill_value=grdROMS.fillval)
                ageice.long_name = "time-averaged age of the ice"
                ageice.units = "years"
                ageice.time = "ocean_time"
                ageice.field = "ice age, scalar, series"
                #ageice.missing_value = grdROMS.fillval

                uice = f1.createVariable('uice', 'd', ('ocean_time', 'eta_u', 'xi_u',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)
                uice.long_name = "time-averaged u-component of ice velocity"
                uice.units = "meter second-1"
                uice.time = "ocean_time"
                uice.field = "u-component of ice velocity, scalar, series"
                #uice.missing_value = grdROMS.fillval

                vice = f1.createVariable('vice', 'd', ('ocean_time', 'eta_v', 'xi_v',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)

                vice.long_name = "time-averaged v-component of ice velocity"
                vice.units = "meter second-1"
                vice.time = "ocean_time"
                vice.field = "v-component of ice velocity, scalar, series"
                #vice.missing_value = grdROMS.fillval

                aice = f1.createVariable('aice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)

                aice.long_name = "time-averaged fraction of cell covered by ice"
                aice.time = "ocean_time"
                aice.field = "ice concentration, scalar, series"
                #aice.missing_value = grdROMS.fillval

                hice = f1.createVariable('hice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)
                hice.long_name = "time-averaged average ice thickness in cell"
                hice.units = "meter"
                hice.time = "ocean_time"
                hice.field = "ice thickness, scalar, series"
                #hice.missing_value = grdROMS.fillval

                snow_thick = f1.createVariable('snow_thick', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                               fill_value=grdROMS.fillval)

                snow_thick.long_name = "time-averaged thickness of snow cover"
                snow_thick.units = "meter"
                snow_thick.time = "ocean_time"
                snow_thick.field = "snow thickness, scalar, series"
                #snow_thick.missing_value = grdROMS.fillval

                ti = f1.createVariable('ti', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)

                ti.long_name = "time-averaged interior ice temperature"
                ti.units = "degrees Celcius"
                ti.time = "ocean_time"
                ti.field = "interior temperature, scalar, series"
                #ti.missing_value = grdROMS.fillval

                sfwat = f1.createVariable('sfwat', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                          fill_value=grdROMS.fillval)

                sfwat.long_name = "time-averaged surface melt water thickness on ice"
                sfwat.units = "meter"
                sfwat.time = "ocean_time"
                sfwat.field = "melt water thickness, scalar, series"
                #sfwat.missing_value = grdROMS.fillval

                tisrf = f1.createVariable('tisrf', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                          fill_value=grdROMS.fillval)
                tisrf.long_name = "time-averaged temperature of ice surface"
                tisrf.units = "degrees Celcius"
                tisrf.time = "ocean_time"
                tisrf.field = "surface temperature, scalar, series"
                #tisrf.missing_value = grdROMS.fillval

                sig11 = f1.createVariable('sig11', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                          fill_value=grdROMS.fillval)

                sig11.long_name = "time-averaged internal ice stress 11 component"
                sig11.units = "Newton meter-1"
                sig11.time = "ocean_time"
                sig11.field = "ice stress 11, scalar, series"
                #sig11.missing_value = grdROMS.fillval

                sig12 = f1.createVariable('sig12', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                          fill_value=grdROMS.fillval)

                sig12.long_name = "time-averaged internal ice stress 12 component"
                sig12.units = "Newton meter-1"
                sig12.time = "ocean_time"
                sig12.field = "ice stress 12, scalar, series"
                #sig12.missing_value = grdROMS.fillval

                sig22 = f1.createVariable('sig22', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                          fill_value=grdROMS.fillval)

                sig22.long_name = "time-averaged internal ice stress 22 component"
                sig22.units = "Newton meter-1"
                sig22.time = "ocean_time"
                sig22.field = "ice stress 22, scalar, series"
                #sig22.missing_value = grdROMS.fillval
            if confM2R.write_bcg:

                v_o3_c = f1.createVariable('O3_c', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_o3_c.long_name = "carbonate/total dissolved inorganic carbon"
                v_o3_c.time = "ocean_time"
                v_o3_c.units = "mmol C/m^3"
                v_o3_c.field = "O3_c, scalar, series"

                v_o3_ta = f1.createVariable('O3_TA', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_o3_ta.long_name = "carbonate/bioalkalinity"
                v_o3_ta.time = "ocean_time"
                v_o3_ta.units = "umol/kg"
                v_o3_ta.field = "O3_ta, scalar, series"

                v_n1_p = f1.createVariable('N1_p', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_n1_p.long_name = "phosphate/phosphorus"
                v_n1_p.time = "ocean_time"
                v_n1_p.units = "mmol P/m^3"
                v_n1_p.field = "N1_p, scalar, series"

                v_o2_o = f1.createVariable('O2_o', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_o2_o.long_name = "oxygen/oxygen"
                v_o2_o.time = "ocean_time"
                v_o2_o.units = "mmol O_2/m^3"
                v_o2_o.field = "O2_o, scalar, series"

                v_n3_n = f1.createVariable('N3_n', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_n3_n.long_name = "nitrate/nitrogen"
                v_n3_n.time = "ocean_time"
                v_n3_n.units = "mmol N/m^3"
                v_n3_n.field = "N3_n, scalar, series"

                v_n5_s = f1.createVariable('N5_s', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
                v_n5_s.long_name = "silicate/silicate"
                v_n5_s.time = "ocean_time"
                v_n5_s.units = "mmol Si/m^3"
                v_n5_s.field = "N5_s, scalar, series"

        # If we are creating climatology files with loops every 360 days, then create these variables here
        if confM2R.isclimatology:
            v_time = f1.createVariable('clim_time', 'd', ('clim_time',), zlib=myzlib, fill_value=grdROMS.fillval)
            v_time.units = 'day'
            v_time.field = 'time, scalar, series'
            v_time.calendar = 'standard'
            v_time.cycle_length = 360.

            v_salt = f1.createVariable('salt', 'f', ('clim_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_salt.long_name = "salinity"
            v_salt.time = "clim_time"
            v_salt.field = "salinity, scalar, series"
            #v_salt.missing_value = grdROMS.fillval

            v_salt = f1.createVariable('SSS', 'f', ('clim_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_salt.long_name = "salinity"
            v_salt.time = "clim_time"
            v_salt.field = "salinity, scalar, series"
            #v_salt.missing_value = grdROMS.fillval

            v_temp = f1.createVariable('temp', 'f', ('clim_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
            v_temp.long_name = "potential temperature"
            v_temp.units = "Celsius"
            v_temp.time = "clim_time"
            v_temp.field = "temperature, scalar, series"
            #v_temp.missing_value = grdROMS.fillval

    else:
        f1 = Dataset(confM2R.clim_name, mode='a', format=confM2R.output_format)

    if confM2R.isclimatology is False:
        if myvar == confM2R.global_varnames[0]:
            print("time units ", grdROMS.timeunits)

            unit = grdROMS.timeunits.split(" ")[0]
            # Convert the units from the datasource into seconds
            secs_per_unit = {"seconds": 1, "hours": 3600, "days": 86400}[unit]
            f1.variables['ocean_time'][ntime] = grdROMS.time * secs_per_unit

            d = num2date(grdROMS.time * secs_per_unit, units=grdROMS.timeunits,
                         calendar=f1.variables['ocean_time'].calendar)
            grdROMS.message = d
            print("Identified secs_per_unit: {} date: {} time: {}".format(secs_per_unit, d, grdROMS.time))

        if myvar == 'temperature':
            f1.variables['temp'][ntime, :, :, :] = data1
        if myvar == 'salinity':
            f1.variables['salt'][ntime, :, :, :] = data1
        if myvar == 'ssh':
            f1.variables['zeta'][ntime, :, :] = data1
        if myvar == 'vvel':
            f1.variables['u'][ntime, :, :, :] = data1
            f1.variables['v'][ntime, :, :, :] = data2

            f1.variables['ubar'][ntime, :, :] = data3
            f1.variables['vbar'][ntime, :, :] = data4

        if confM2R.write_ice:
            if myvar == "ageice":
                # print "NOTE! Setting values of ageice to ZERO! (IOWrite.py)"
                data1 = np.where(abs(data1) > 100, 0, data1)
                f1.variables['ageice'][ntime, :, :] = data1

            if myvar == 'uice':
                data1 = np.where(abs(data1) > 120, 0, data1)
                f1.variables['uice'][ntime, :, :] = data1 * 0.01  # NorESM is cm/s divide by 100 to get m/s
                f1.variables['sfwat'][ntime, :, :] = 0.
                f1.variables['tisrf'][ntime, :, :] = 0.
                f1.variables['ti'][ntime, :, :] = 0.
                f1.variables['sig11'][ntime, :, :] = 0.
                f1.variables['sig12'][ntime, :, :] = 0.
                f1.variables['sig22'][ntime, :, :] = 0.

                if confM2R.ocean_indata_type == 'GLORYS':
                    # Special care for GLORYS as dataset does not contain sea ice age and snow thickness
                    f1.variables['ageice'][ntime, :, :] = 0.
                    f1.variables['snow_thick'][ntime, :, :] = 0

            if myvar == 'vice':
                data1 = np.where(abs(data1) > 120, 0, data1)
                f1.variables['vice'][ntime, :, :] = data1 * 0.01  # NorESM is cm/s divide by 100 to get m/s
            if myvar == 'aice':
                data1 = np.where(abs(data1) > 120, 0, data1)
                f1.variables['aice'][ntime, :, :] = data1 * 0.01  # NorESM is % divide by 100 to get fraction
            if myvar == 'hice':
                data1 = np.where(abs(data1) > 10, 0, data1)
                # data1 = np.ma.masked_where(abs(data1) > 10, data1)
                f1.variables['hice'][ntime, :, :] = data1
            if myvar == 'snow_thick':
                # data1 = np.ma.masked_where(abs(data1) > 100, data1)
                data1 = np.where(abs(data1) > 10, 0, data1)
                f1.variables['snow_thick'][ntime, :, :] = data1

        if confM2R.write_bcg:
            if myvar in ['O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
                data1 = np.where(abs(data1) < 0, 0, data1)
                if confM2R.ocean_indata_type== "NORESM":
                    """
                    Multiply the NORESM variable by conversion factors below:
                    NORESM name     NORESM units     ERSEM name     ERSEM units    Conversion factor
                    dissic                      [mol C/m3]           O3_c                   [mmol C/m3]   1e3
                    talk                         [eq/m3]                 O3_TA                [umol/kg]         1e6/1025
                    po4                         [mol P/m3]           N1_p                   [mmol P/m3]   1e3
                    no3                         [mol N/m3]           N3_n                  [mmol N/m3]   1e3
                    si                             [mol Si/m3]           N5_s                   [mmol Si/m3]   1e3
                    o2                           [mol O2/m3]         O2_o                  [mmol O2/m3] 1e3
                    """
                
                    if myvar=="O3_TA": 
                        data1=data1*1.0e6/1025.
                    else: 
                        data1=data1*1.0e3
                f1.variables[myvar][ntime,:,:,:] = data1

    if confM2R.isclimatology:
        # Climatological time starts at the 15th of each month
        d = datetime(2012, int(ntime) + 1, 1)
        tt = d.timetuple()
        if myvar == "temperature":
            f1.variables['clim_time'][ntime] = tt.tm_yday + 15

        grdROMS.message = tt.tm_yday + 15

        if myvar == 'temperature':
            f1.variables['temp'][ntime,:,:,:] = data1
        if myvar == 'salinity':
            f1.variables['salt'][ntime,:,:,:] = data1
            f1.variables['SSS'][ntime,:,:,:] = data1

    f1.close()
