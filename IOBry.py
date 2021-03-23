import os
import time
from datetime import datetime

from netCDF4 import Dataset

_author_ = 'Trond Kristiansen'
_email_ = 'me@trondkristiansen.com'
_created_ = datetime(2009, 3, 2)
_modified_ = datetime(2014, 4, 7)
_version_ = "0.1.0"
_status_ = "Development"


def help():
    """
    This function generates a BRY file from scratch. The variables are all created
    for East, West, North, and South. Varibales include:
    salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is ocean_time which is days
    since 1948/1/1.

    This file is netcdf CF compliant and follows the setup for variable names and units given in the ROMS source
    file: Data/ROMS/CDL/bry_unlimit.cdl

    (also see: https://www.myroms.org/forum/viewtopic.php?f=30&t=1450&p=5209&hilit=cf+compliant#p5209)

    This function is called from clim2bry.py.

    To check the BRY file for CF compliancy: http://titania.badc.rl.ac.uk/cgi-bin/cf-checker.pl?cfversion=1.0
    """


def createBryFile(confM2R):
    if (confM2R.output_format == 'NETCDF4'):
        myzlib = True
    else:
        myzlib = False
    grdROMS = confM2R.grdROMS

    if os.path.exists(confM2R.bry_name):
        os.remove(confM2R.bry_name)
    print(('\n=>Creating initial Boundary (BRY) file {}'.format(confM2R.bry_name)))

    f1 = Dataset(confM2R.bry_name, mode='w', format=confM2R.output_format)
    f1.title = "Boundary forcing file (BRY) used for forcing of the ROMS model"
    f1.description = "Created for the {} grid file".format(confM2R.roms_grid_path)
    f1.grdFile = "{}".format(confM2R.roms_grid_path)
    f1.history = 'Created ' + time.ctime(time.time())
    f1.source = "{} ({})".format(confM2R.authorname, confM2R.authoremail)
    f1.type = "File in {} format created using MODEL2ROMS".format(confM2R.output_format)
    f1.link = "https://github.com/trondkr/model2roms"
    f1.Conventions = "CF-1.0"

    """ Define dimensions """
    f1.createDimension('xi_rho', grdROMS.xi_rho)
    f1.createDimension('eta_rho', grdROMS.eta_rho)
    f1.createDimension('xi_u', grdROMS.xi_u)
    f1.createDimension('eta_u', grdROMS.eta_u)
    f1.createDimension('xi_v', grdROMS.xi_v)
    f1.createDimension('eta_v', grdROMS.eta_v)
    f1.createDimension('xi_psi', grdROMS.xi_psi)
    f1.createDimension('eta_psi', grdROMS.eta_psi)
    f1.createDimension('ocean_time', None)
    f1.createDimension('s_rho', len(grdROMS.s_rho))
    f1.createDimension('s_w', len(grdROMS.s_w))

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
    vnc.long_name = "Bathymetry at RHO-points"
    vnc.units = "meter"
    vnc.coordinates = "lat_rho lon_rho"
    vnc.field = "bath, scalar"
    vnc[:, :] = grdROMS.h

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
    vnc.field = "Cs_rho, scalar"
    vnc[:] = grdROMS.Cs_rho

    vnc = f1.createVariable('Cs_w', 'd', ('s_w',), zlib=myzlib, fill_value=grdROMS.fillval)
    vnc.long_name = "S-coordinate stretching curves at W-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "Cs_w, scalar"
    vnc[:] = grdROMS.Cs_w

    vnc = f1.createVariable('hc', 'd')
    vnc.long_name = "S-coordinate parameter, critical depth";
    vnc.units = "meter"
    vnc[:] = grdROMS.hc

    vnc = f1.createVariable('z_r', 'd', ('s_rho', 'eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
    vnc.long_name = "Sigma layer to depth matrix";
    vnc.units = "meter"
    vnc[:, :, :] = grdROMS.z_r

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

    v_time = f1.createVariable('ocean_time', 'd', ('ocean_time',), zlib=myzlib, fill_value=grdROMS.fillval)
    
    v_time.long_name = 'seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'

    if (confM2R.ocean_indata_type == "NORESM"):
        v_time.calendar = 'noleap'
    else:
        v_time.calendar = 'standard'

    v_temp_west = f1.createVariable('temp_west', 'f', ('ocean_time', 's_rho', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_temp_west.long_name = "potential temperature western boundary condition"
    v_temp_west.units = "Celsius"
    v_temp_west.field = "temp_west, scalar, series"
    #v_temp_west.missing_value = grdROMS.fillval
    v_temp_west.time = "ocean_time"

    v_temp_east = f1.createVariable('temp_east', 'f', ('ocean_time', 's_rho', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_temp_east.long_name = "potential temperature eastern boundary condition"
    v_temp_east.units = "Celsius"
    v_temp_east.field = "temp_east, scalar, series"
    #v_temp_east.missing_value = grdROMS.fillval
    v_temp_east.time = "ocean_time"

    v_temp_south = f1.createVariable('temp_south', 'f', ('ocean_time', 's_rho', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_temp_south.long_name = "potential temperature southern boundary condition"
    v_temp_south.units = "Celsius"
    v_temp_south.field = "temp_south, scalar, series"
    #v_temp_south.missing_value = grdROMS.fillval
    v_temp_south.time = "ocean_time"

    v_temp_north = f1.createVariable('temp_north', 'f', ('ocean_time', 's_rho', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_temp_north.long_name = "potential temperature northern boundary condition"
    v_temp_north.units = "Celsius"
    v_temp_north.field = "temp_north, scalar, series"
    #v_temp_north.missing_value = grdROMS.fillval
    v_temp_north.time = "ocean_time"

    v_salt_west = f1.createVariable('salt_west', 'f', ('ocean_time', 's_rho', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_salt_west.long_name = "salinity western boundary condition"
    v_salt_west.field = "salt_west, scalar, series"
    #v_salt_west.missing_value = grdROMS.fillval
    v_salt_west.time = "ocean_time"

    v_salt_east = f1.createVariable('salt_east', 'f', ('ocean_time', 's_rho', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_salt_east.long_name = "salinity eastern boundary condition"
    v_salt_east.field = "salt_east, scalar, series"
    #v_salt_east.missing_value = grdROMS.fillval
    v_salt_east.time = "ocean_time"

    v_salt_south = f1.createVariable('salt_south', 'f', ('ocean_time', 's_rho', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_salt_south.long_name = "salinity southern boundary condition"
    v_salt_south.field = "salt_south, scalar, series"
    #v_salt_south.missing_value = grdROMS.fillval
    v_salt_south.time = "ocean_time"

    v_salt_north = f1.createVariable('salt_north', 'f', ('ocean_time', 's_rho', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_salt_north.long_name = "salinity northern boundary condition"
    v_salt_north.field = "salt_north, scalar, series"
    #v_salt_north.missing_value = grdROMS.fillval
    v_salt_north.time = "ocean_time"

    v_ssh_west = f1.createVariable('zeta_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                   fill_value=grdROMS.fillval)
    v_ssh_west.long_name = "free-surface western boundary condition"
    v_ssh_west.units = "meter"
    v_ssh_west.field = "zeta_west, scalar, series"
    #v_ssh_west.missing_value = grdROMS.fillval
    v_ssh_west.time = "ocean_time"

    v_ssh_east = f1.createVariable('zeta_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                   fill_value=grdROMS.fillval)
    v_ssh_east.long_name = "free-surface eastern boundary condition"
    v_ssh_east.units = "meter"
    v_ssh_east.field = "zeta_east, scalar, series"
    #v_ssh_east.missing_value = grdROMS.fillval
    v_ssh_east.time = "ocean_time"

    v_ssh_south = f1.createVariable('zeta_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_ssh_south.long_name = "free-surface southern boundary condition"
    v_ssh_south.units = "meter"
    v_ssh_south.field = "zeta_south, scalar, series"
    #v_ssh_south.missing_value = grdROMS.fillval
    v_ssh_south.time = "ocean_time"

    v_ssh_north = f1.createVariable('zeta_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_ssh_north.long_name = "free-surface northern boundary condition"
    v_ssh_north.units = "meter"
    v_ssh_north.field = "zeta_north, scalar, series"
    #v_ssh_north.missing_value = grdROMS.fillval
    v_ssh_north.time = "ocean_time"

    v_u_west = f1.createVariable('u_west', 'f', ('ocean_time', 's_rho', 'eta_u',), zlib=myzlib,
                                 fill_value=grdROMS.fillval)
    v_u_west.long_name = "3D u-momentum western boundary condition"
    v_u_west.units = "meter second-1"
    v_u_west.field = "u_west, scalar, series"
    #v_u_west.missing_value = grdROMS.fillval
    v_u_west.time = "ocean_time"

    v_u_east = f1.createVariable('u_east', 'f', ('ocean_time', 's_rho', 'eta_u',), zlib=myzlib,
                                 fill_value=grdROMS.fillval)
    v_u_east.long_name = "3D u-momentum eastern boundary condition"
    v_u_east.units = "meter second-1"
    v_u_east.field = "u_east, scalar, series"
    #v_u_east.missing_value = grdROMS.fillval
    v_u_east.time = "ocean_time"

    v_u_south = f1.createVariable('u_south', 'f', ('ocean_time', 's_rho', 'xi_u',), zlib=myzlib,
                                  fill_value=grdROMS.fillval)
    v_u_south.long_name = "3D u-momentum southern boundary condition"
    v_u_south.units = "meter second-1"
    v_u_south.field = "u_south, scalar, series"
    #v_u_south.missing_value = grdROMS.fillval
    v_u_south.time = "ocean_time"

    v_u_north = f1.createVariable('u_north', 'f', ('ocean_time', 's_rho', 'xi_u',), zlib=myzlib,
                                  fill_value=grdROMS.fillval)
    v_u_north.long_name = "3D u-momentum northern boundary condition"
    v_u_north.units = "meter second-1"
    v_u_north.field = "u_north, scalar, series"
    #v_u_north.missing_value = grdROMS.fillval
    v_u_north.time = "ocean_time"

    v_v_west = f1.createVariable('v_west', 'f', ('ocean_time', 's_rho', 'eta_v',), zlib=myzlib,
                                 fill_value=grdROMS.fillval)
    v_v_west.long_name = "3D v-momentum western boundary condition"
    v_v_west.units = "meter second-1"
    v_v_west.field = "v_west, scalar, series"
    #v_v_west.missing_value = grdROMS.fillval
    v_v_west.time = "ocean_time"

    v_v_east = f1.createVariable('v_east', 'f', ('ocean_time', 's_rho', 'eta_v',), zlib=myzlib,
                                 fill_value=grdROMS.fillval)
    v_v_east.long_name = "3D v-momentum eastern boundary condition"
    v_v_east.units = "meter second-1"
    v_v_east.field = "v_east, scalar, series"
    #v_v_east.missing_value = grdROMS.fillval
    v_v_east.time = "ocean_time"

    v_v_south = f1.createVariable('v_south', 'f', ('ocean_time', 's_rho', 'xi_v',), zlib=myzlib,
                                  fill_value=grdROMS.fillval)
    v_v_south.long_name = "3D v-momentum southern boundary condition"
    v_v_south.units = "meter second-1"
    v_v_south.field = "v_south, scalar, series"
    #v_v_south.missing_value = grdROMS.fillval
    v_v_south.time = "ocean_time"

    v_v_north = f1.createVariable('v_north', 'f', ('ocean_time', 's_rho', 'xi_v',), zlib=myzlib,
                                  fill_value=grdROMS.fillval)
    v_v_north.long_name = "3D v-momentum northern boundary condition"
    v_v_north.units = "meter second-1"
    v_v_north.field = "v_north, scalar, series"
    #v_v_north.missing_value = grdROMS.fillval
    v_v_north.time = "ocean_time"

    v_vbar_west = f1.createVariable('vbar_west', 'f', ('ocean_time', 'eta_v',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_vbar_west.long_name = "2D v-momentum western boundary condition"
    v_vbar_west.units = "meter second-1"
    v_vbar_west.field = "vbar_west, scalar, series"
    #v_vbar_west.missing_value = grdROMS.fillval
    v_vbar_west.time = "ocean_time"

    v_vbar_east = f1.createVariable('vbar_east', 'f', ('ocean_time', 'eta_v',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_vbar_east.long_name = "2D v-momentum eastern boundary condition"
    v_vbar_east.units = "meter second-1"
    v_vbar_east.field = "vbar_east, scalar, series"
    #v_vbar_east.missing_value = grdROMS.fillval
    v_vbar_east.time = "ocean_time"

    v_vbar_south = f1.createVariable('vbar_south', 'f', ('ocean_time', 'xi_v',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_vbar_south.long_name = "2D v-momentum southern boundary condition"
    v_vbar_south.units = "meter second-1"
    v_vbar_south.field = "vbar_south, scalar, series"
    #v_vbar_south.missing_value = grdROMS.fillval
    v_vbar_south.time = "ocean_time"

    v_vbar_north = f1.createVariable('vbar_north', 'f', ('ocean_time', 'xi_v',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_vbar_north.long_name = "2D v-momentum northern boundary condition"
    v_vbar_north.units = "meter second-1"
    v_vbar_north.field = "vbar_north, scalar, series"
    #v_vbar_north.missing_value = grdROMS.fillval
    v_vbar_north.time = "ocean_time"

    v_ubar_west = f1.createVariable('ubar_west', 'f', ('ocean_time', 'eta_u',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_ubar_west.long_name = "2D u-momentum western boundary condition"
    v_ubar_west.units = "meter second-1"
    v_ubar_west.field = "ubar_west, scalar, series"
   # v_ubar_west.missing_value = grdROMS.fillval
    v_ubar_west.time = "ocean_time"

    v_ubar_east = f1.createVariable('ubar_east', 'f', ('ocean_time', 'eta_u',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
    v_ubar_east.long_name = "2D u-momentum eastern boundary condition"
    v_ubar_east.units = "meter second-1"
    v_ubar_east.field = "ubar_east, scalar, series"
    #v_ubar_east.missing_value = grdROMS.fillval
    v_ubar_east.time = "ocean_time"

    v_ubar_south = f1.createVariable('ubar_south', 'f', ('ocean_time', 'xi_u',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_ubar_south.long_name = "2D u-momentum southern boundary condition"
    v_ubar_south.units = "meter second-1"
    v_ubar_south.field = "ubar_south, scalar, series"
    #v_ubar_south.missing_value = grdROMS.fillval
    v_ubar_south.time = "ocean_time"

    v_ubar_north = f1.createVariable('ubar_north', 'f', ('ocean_time', 'xi_u',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
    v_ubar_north.long_name = "2D u-momentum northern boundary condition"
    v_ubar_north.units = "meter second-1"
    v_ubar_north.field = "ubar_north, scalar, series"
    #v_ubar_north.missing_value = grdROMS.fillval
    v_ubar_north.time = "ocean_time"

    if confM2R.write_bcg:
        directions=['east','west','north','south']
        dimens=['eta_rho','eta_rho','xi_rho','xi_rho']

        lndirections=['eastern','western','northern','southern']
        for currentdir,lndirection,dim in zip(directions, lndirections,dimens):
            currentvar='O3_c_'+currentdir
            longname="carbonate/total dissolved inorganic carbon {} boundary condition".format(lndirection) 
            O3_c = f1.createVariable(currentvar, 'f', ('ocean_time', 's_rho', dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            O3_c.long_name = longname
            O3_c.field = "{}, scalar, series".format(currentvar) 
            O3_c.units = "mmol C/m^3" 

            currentvar='O3_TA_'+currentdir
            longname="carbonate/bioalkalinity {} boundary condition".format(lndirection) 
            O3_ta = f1.createVariable(currentvar, 'f', ('ocean_time', 's_rho', dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            O3_ta.long_name = longname
            O3_ta.field = "{}, scalar, series".format(currentvar) 
            O3_ta.units = "umol/kg" 

            currentvar='N1_p_'+currentdir
            longname="phosphate/phosphorus {} boundary condition".format(lndirection) 
            N1_p = f1.createVariable(currentvar, 'f', ('ocean_time','s_rho', dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            N1_p.long_name = longname
            N1_p.field = "{}, scalar, series".format(currentvar) 
            N1_p.units = "mmol P/m^3" 

            currentvar='N3_n_'+currentdir
            longname="nitrate/nitrogen {} boundary condition".format(lndirection) 
            N3_n = f1.createVariable(currentvar, 'f', ('ocean_time','s_rho', dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            N3_n.long_name = longname
            N3_n.field = "{}, scalar, series".format(currentvar) 
            N3_n.units = "mmol N/m^3" 

            currentvar='N5_s_'+currentdir
            longname="silicate/silicate {} boundary condition".format(lndirection) 
            N5_s = f1.createVariable(currentvar, 'f', ('ocean_time','s_rho', dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            N5_s.long_name = longname
            N5_s.field = "{}, scalar, series".format(currentvar) 
            N5_s.units = "mmol Si/m^3" 

            currentvar='O2_o_'+currentdir
            longname="oxygen/oxygen {} boundary condition".format(lndirection) 
            O2_o = f1.createVariable(currentvar, 'f', ('ocean_time', 's_rho',dim,), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
            O2_o.long_name = longname
            O2_o.field = "{}, scalar, series".format(currentvar) 
            O2_o.units = "mmol O_2/m^3" 

    if confM2R.write_ice:
        ageice_west = f1.createVariable('ageice_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        ageice_west.long_name = "time-averaged age of the ice western boundary conditions"
        ageice_west.units = "years"
        ageice_west.time = "ocean_time"
        ageice_west.field = "ice age, scalar, series"
        #ageice_west.missing_value = grdROMS.fillval

        ageice_east = f1.createVariable('ageice_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        ageice_east.long_name = "time-averaged age of the ice eastern boundary conditions"
        ageice_east.units = "years"
        ageice_east.time = "ocean_time"
        ageice_east.field = "ice age, scalar, series"
        #ageice_east.missing_value = grdROMS.fillval

        ageice_south = f1.createVariable('ageice_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)
        ageice_south.long_name = "time-averaged age of the ice southern boundary conditions"
        ageice_south.units = "years"
        ageice_south.time = "ocean_time"
        ageice_south.field = "ice age, scalar, series"
        #ageice_south.missing_value = grdROMS.fillval

        ageice_north = f1.createVariable('ageice_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                         fill_value=grdROMS.fillval)
        ageice_north.long_name = "time-averaged age of the ice northern boundary conditions"
        ageice_north.units = "years"
        ageice_north.time = "ocean_time"
        ageice_north.field = "ice age, scalar, series"
        #ageice_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        uice_west = f1.createVariable('uice_west', 'f', ('ocean_time', 'eta_u',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        uice_west.long_name = "time-averaged age of the u-component of ice velocity western boundary conditions"
        uice_west.units = "meter second-1"
        uice_west.time = "ocean_time"
        uice_west.field = "u-component of ice velocity, scalar, series"
        #uice_west.missing_value = grdROMS.fillval

        uice_east = f1.createVariable('uice_east', 'f', ('ocean_time', 'eta_u',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        uice_east.long_name = "time-averaged age of the u-component of ice velocity eastern boundary conditions"
        uice_east.units = "meter second-1"
        uice_east.time = "ocean_time"
        uice_east.field = "u-component of ice velocity, scalar, series"
        #uice_east.missing_value = grdROMS.fillval

        uice_south = f1.createVariable('uice_south', 'f', ('ocean_time', 'xi_u',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        uice_south.long_name = "time-averaged age of the u-component of ice velocity southern boundary conditions"
        uice_south.units = "meter second-1"
        uice_south.time = "ocean_time"
        uice_south.field = "u-component of ice velocity, scalar, series"
        #uice_south.missing_value = grdROMS.fillval

        uice_north = f1.createVariable('uice_north', 'f', ('ocean_time', 'xi_u',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        uice_north.long_name = "time-averaged age of the u-component of ice velocity northern boundary conditions"
        uice_north.units = "meter second-1"
        uice_north.time = "ocean_time"
        uice_north.field = "u-component of ice velocity, scalar, series"
        #uice_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        vice_west = f1.createVariable('vice_west', 'f', ('ocean_time', 'eta_v',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        vice_west.long_name = "time-averaged age of the v-component of ice velocity western boundary conditions"
        vice_west.units = "meter second-1"
        uice_west.time = "ocean_time"
        vice_west.field = "v-component of ice velocity, scalar, series"
        #vice_west.missing_value = grdROMS.fillval

        vice_east = f1.createVariable('vice_east', 'f', ('ocean_time', 'eta_v',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        vice_east.long_name = "time-averaged age of the v-component of ice velocity eastern boundary conditions"
        vice_east.units = "meter second-1"
        vice_east.time = "ocean_time"
        vice_east.field = "v-component of ice velocity, scalar, series"
        #vice_east.missing_value = grdROMS.fillval

        vice_south = f1.createVariable('vice_south', 'f', ('ocean_time', 'xi_v',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        vice_south.long_name = "time-averaged age of the v-component of ice velocity southern boundary conditions"
        vice_south.units = "meter second-1"
        vice_south.time = "ocean_time"
        vice_south.field = "v-component of ice velocity, scalar, series"
        #vice_south.missing_value = grdROMS.fillval

        vice_north = f1.createVariable('vice_north', 'f', ('ocean_time', 'xi_v',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        vice_north.long_name = "time-averaged age of the u-component of ice velocity northern boundary conditions"
        vice_north.units = "meter second-1"
        vice_north.time = "ocean_time"
        vice_north.field = "v-component of ice velocity, scalar, series"
        #vice_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        aice_west = f1.createVariable('aice_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        aice_west.long_name = "time-averaged fraction of cell covered by ice western boundary conditions"
        aice_west.units = "%"
        aice_west.time = "ocean_time"
        aice_west.field = "ice concentration, scalar, series"
        #aice_west.missing_value = grdROMS.fillval

        aice_east = f1.createVariable('aice_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        aice_east.long_name = "time-averaged fraction of cell covered by ice eastern boundary conditions"
        aice_east.units = "%"
        aice_east.time = "ocean_time"
        aice_east.field = "ice concentration, scalar, series"
        #aice_east.missing_value = grdROMS.fillval

        aice_south = f1.createVariable('aice_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        aice_south.long_name = "time-averaged fraction of cell covered by ice southern boundary conditions"
        aice_south.units = "%"
        aice_south.time = "ocean_time"
        aice_south.field = "ice concentration, scalar, series"
        #aice_south.missing_value = grdROMS.fillval

        aice_north = f1.createVariable('aice_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        aice_north.long_name = "time-averaged fraction of cell covered by ice northern boundary conditions"
        aice_north.units = "%"
        aice_north.time = "ocean_time"
        aice_north.field = "ice concentration, scalar, series"
        #aice_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        hice_west = f1.createVariable('hice_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        hice_west.long_name = "time-averaged ice thickness in cell western boundary conditions"
        hice_west.units = "meter"
        hice_west.time = "ocean_time"
        hice_west.field = "ice thickness, scalar, series"
        #hice_west.missing_value = grdROMS.fillval

        hice_east = f1.createVariable('hice_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                      fill_value=grdROMS.fillval)
        hice_east.long_name = "time-averaged ice thickness in cell eastern boundary conditions"
        hice_east.units = "meter"
        hice_east.time = "ocean_time"
        hice_east.field = "ice thickness, scalar, series"
        #hice_east.missing_value = grdROMS.fillval

        hice_south = f1.createVariable('hice_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        hice_south.long_name = "time-averaged ice thickness in cell southern boundary conditions"
        hice_south.units = "meter"
        hice_south.time = "ocean_time"
        hice_south.field = "ice thickness, scalar, series"
        #hice_south.missing_value = grdROMS.fillval

        hice_north = f1.createVariable('hice_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        hice_north.long_name = "time-averaged ice thickness in cell northern boundary conditions"
        hice_north.units = "meter"
        hice_north.time = "ocean_time"
        hice_north.field = "ice thickness, scalar, series"
        #hice_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        snow_thick_west = f1.createVariable('snow_thick_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                            fill_value=grdROMS.fillval)
        snow_thick_west.long_name = "time-averaged ice thickness in cell western boundary conditions"
        snow_thick_west.units = "meter"
        snow_thick_west.time = "ocean_time"
        snow_thick_west.field = "snow thickness, scalar, series"
        #snow_thick_west.missing_value = grdROMS.fillval

        snow_thick_east = f1.createVariable('snow_thick_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                            fill_value=grdROMS.fillval)
        snow_thick_east.long_name = "time-averaged ice thickness in cell eastern boundary conditions"
        snow_thick_east.units = "meter"
        snow_thick_east.time = "ocean_time"
        snow_thick_east.field = "snow thickness, scalar, series"
        #snow_thick_east.missing_value = grdROMS.fillval

        snow_thick_south = f1.createVariable('snow_thick_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                             fill_value=grdROMS.fillval)
        snow_thick_south.long_name = "time-averaged ice thickness in cell southern boundary conditions"
        snow_thick_south.units = "meter"
        snow_thick_south.time = "ocean_time"
        snow_thick_south.field = "snow thickness, scalar, series"
        #snow_thick_south.missing_value = grdROMS.fillval

        snow_thick_north = f1.createVariable('snow_thick_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                             fill_value=grdROMS.fillval)
        snow_thick_north.long_name = "time-averaged ice thickness in cell northern boundary conditions"
        snow_thick_north.units = "meter"
        snow_thick_north.time = "ocean_time"
        snow_thick_north.field = "snow thickness, scalar, series"
        #snow_thick_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        ti_west = f1.createVariable('ti_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
        ti_west.long_name = "time-averaged interior ice temperature cell western boundary conditions"
        ti_west.units = "degrees Celcius"
        ti_west.time = "ocean_time"
        ti_west.field = "interior temperature, scalar, series"
        #ti_west.missing_value = grdROMS.fillval

        ti_east = f1.createVariable('ti_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                    fill_value=grdROMS.fillval)
        ti_east.long_name = "time-averaged interior ice temperature eastern boundary conditions"
        ti_east.units = "degrees Celcius"
        ti_east.time = "ocean_time"
        ti_east.field = "interior temperature, scalar, series"
        #ti_east.missing_value = grdROMS.fillval

        ti_south = f1.createVariable('ti_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
        ti_south.long_name = "time-averaged interior ice temperature southern boundary conditions"
        ti_south.units = "degrees Celcius"
        ti_south.time = "ocean_time"
        ti_south.field = "interior temperature, scalar, series"
        #ti_south.missing_value = grdROMS.fillval

        ti_north = f1.createVariable('ti_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                     fill_value=grdROMS.fillval)
        ti_north.long_name = "time-averaged interior ice temperature northern boundary conditions"
        ti_north.units = "degrees Celcius"
        ti_north.time = "ocean_time"
        ti_north.field = "interior temperature, scalar, series"
        #ti_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        sfwat_west = f1.createVariable('sfwat_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sfwat_west.long_name = "time-averaged surface melt water thickness on ice western boundary conditions"
        sfwat_west.units = "meter"
        sfwat_west.time = "ocean_time"
        sfwat_west.field = "melt water thickness, scalar, series"
        #sfwat_west.missing_value = grdROMS.fillval

        sfwat_east = f1.createVariable('sfwat_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sfwat_east.long_name = "time-averaged surface melt water thickness on ice eastern boundary conditions"
        sfwat_east.units = "meter"
        sfwat_east.time = "ocean_time"
        sfwat_east.field = "melt water thickness, scalar, series"
        #sfwat_east.missing_value = grdROMS.fillval

        sfwat_south = f1.createVariable('sfwat_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sfwat_south.long_name = "time-averaged surface melt water thickness on ice southern boundary conditions"
        sfwat_south.units = "meter"
        sfwat_south.time = "ocean_time"
        sfwat_south.field = "melt water thickness, scalar, series"
        #sfwat_south.missing_value = grdROMS.fillval

        sfwat_north = f1.createVariable('sfwat_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sfwat_north.long_name = "time-averaged surface melt water thickness on ice northern boundary conditions"
        sfwat_north.units = "meter"
        sfwat_north.time = "ocean_time"
        sfwat_north.field = "melt water thickness, scalar, series"
        #sfwat_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        tisrf_west = f1.createVariable('tisrf_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        tisrf_west.long_name = "time-averaged temperature of ice surfacewestern boundary conditions"
        tisrf_west.units = "degrees Celcius"
        tisrf_west.time = "ocean_time"
        tisrf_west.field = "surface temperature, scalar, series"
        #tisrf_west.missing_value = grdROMS.fillval

        tisrf_east = f1.createVariable('tisrf_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        tisrf_east.long_name = "time-averaged temperature of ice surface eastern boundary conditions"
        tisrf_east.units = "degrees Celcius"
        tisrf_east.time = "ocean_time"
        tisrf_east.field = "surface temperature, scalar, series"
        #tisrf_east.missing_value = grdROMS.fillval

        tisrf_south = f1.createVariable('tisrf_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        tisrf_south.long_name = "time-averaged temperature of ice surface southern boundary conditions"
        tisrf_south.units = "degrees Celcius"
        tisrf_south.time = "ocean_time"
        tisrf_south.field = "surface temperature, scalar, series"
        #tisrf_south.missing_value = grdROMS.fillval

        tisrf_north = f1.createVariable('tisrf_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        tisrf_north.long_name = "time-averaged temperature of ice surface northern boundary conditions"
        tisrf_north.units = "degrees Celcius"
        tisrf_north.time = "ocean_time"
        tisrf_north.field = "surface temperature, scalar, series"
        #tisrf_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        sig11_west = f1.createVariable('sig11_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig11_west.long_name = "time-averaged internal ice stress 11 component boundary conditions"
        sig11_west.units = "Newton meter-1"
        sig11_west.time = "ocean_time"
        sig11_west.field = "ice stress 11, scalar, series"
        #sig11_west.missing_value = grdROMS.fillval

        sig11_east = f1.createVariable('sig11_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig11_east.long_name = "time-averaged internal ice stress 11 component eastern boundary conditions"
        sig11_east.units = "Newton meter-1"
        sig11_east.time = "ocean_time"
        sig11_east.field = "ice stress 11, scalar, series"
        #sig11_east.missing_value = grdROMS.fillval

        sig11_south = f1.createVariable('sig11_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig11_south.long_name = "time-averaged internal ice stress 11 componentsouthern boundary conditions"
        sig11_south.units = "Newton meter-1"
        sig11_south.time = "ocean_time"
        sig11_south.field = "ice stress 11, scalar, series"
        #sig11_south.missing_value = grdROMS.fillval

        sig11_north = f1.createVariable('sig11_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig11_north.long_name = "time-averaged internal ice stress 11 component northern boundary conditions"
        sig11_north.units = "Newton meter-1"
        sig11_north.time = "ocean_time"
        sig11_north.field = "ice stress 11, scalar, series"
        #sig11_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        sig12_west = f1.createVariable('sig12_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig12_west.long_name = "time-averaged internal ice stress 12 component boundary conditions"
        sig12_west.units = "Newton meter-1"
        sig12_west.time = "ocean_time"
        sig12_west.field = "ice stress 12, scalar, series"
        #sig12_west.missing_value = grdROMS.fillval

        sig12_east = f1.createVariable('sig12_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig12_east.long_name = "time-averaged internal ice stress 12 component eastern boundary conditions"
        sig12_east.units = "Newton meter-1"
        sig12_east.time = "ocean_time"
        sig12_east.field = "ice stress 12, scalar, series"
        #sig12_east.missing_value = grdROMS.fillval

        sig12_south = f1.createVariable('sig12_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig12_south.long_name = "time-averaged internal ice stress 12 componentsouthern boundary conditions"
        sig12_south.units = "Newton meter-1"
        sig12_south.time = "ocean_time"
        sig12_south.field = "ice stress 12, scalar, series"
        #sig12_south.missing_value = grdROMS.fillval

        sig12_north = f1.createVariable('sig12_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig12_north.long_name = "time-averaged internal ice stress 12 component northern boundary conditions"
        sig12_north.units = "Newton meter-1"
        sig12_north.time = "ocean_time"
        sig12_north.field = "ice stress 12, scalar, series"
        #sig12_north.missing_value = grdROMS.fillval

        # ----------------------------------------

        sig22_west = f1.createVariable('sig22_west', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig22_west.long_name = "time-averaged internal ice stress 22 component boundary conditions"
        sig22_west.units = "Newton meter-1"
        sig22_west.time = "ocean_time"
        sig22_west.field = "ice stress 22, scalar, series"
        #sig22_west.missing_value = grdROMS.fillval

        sig22_east = f1.createVariable('sig22_east', 'f', ('ocean_time', 'eta_rho',), zlib=myzlib,
                                       fill_value=grdROMS.fillval)
        sig22_east.long_name = "time-averaged internal ice stress 22 component eastern boundary conditions"
        sig22_east.units = "Newton meter-1"
        sig22_east.time = "ocean_time"
        sig22_east.field = "ice stress 22, scalar, series"
        #sig22_east.missing_value = grdROMS.fillval

        sig22_south = f1.createVariable('sig22_south', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig22_south.long_name = "time-averaged internal ice stress 22 componentsouthern boundary conditions"
        sig22_south.units = "Newton meter-1"
        sig22_south.time = "ocean_time"
        sig22_south.field = "ice stress 22, scalar, series"
        #sig22_south.missing_value = grdROMS.fillval

        sig22_north = f1.createVariable('sig22_north', 'f', ('ocean_time', 'xi_rho',), zlib=myzlib,
                                        fill_value=grdROMS.fillval)
        sig22_north.long_name = "time-averaged internal ice stress 22 component northern boundary conditions"
        sig22_north.units = "Newton meter-1"
        sig22_north.time = "ocean_time"
        sig22_north.field = "ice stress 22, scalar, series"
        #sig22_north.missing_value = grdROMS.fillval

        # ----------------------------------------

    f1.close()
