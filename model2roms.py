from __future__ import print_function

import logging
from datetime import datetime

import barotropic
import interpolation as interp
import numpy as np
from netCDF4 import Dataset, date2num

import IOinitial
import IOsubset
import IOwrite
import datetimeFunctions
import forcingFilenames as fc
import interp2D

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass
__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2008, 8, 15)
__modified__ = datetime(2021, 3, 23)
__version__ = "1.8"
__status__ = "Development, modified on 15.08.2008,01.10.2009,07.01.2010, " \
             "15.07.2014, 01.12.2014, 07.08.2015, " \
             "08.02.2018, 04.03.2019, 13.03.2019, 23.03.2021"


def vertical_interpolation(myvar, array1, array2, grdROMS, grdMODEL):
    outINDEX_ST = (grdROMS.nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    outINDEX_U = (grdROMS.nlevels, grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_V = (grdROMS.nlevels, grdROMS.eta_v, grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v, grdROMS.xi_v)

    if myvar in ['salinity', 'temperature', 'O3_c', 'O3_TA', 'N1_p', 'N3_n', 'N5_s', 'O2_o']:
        logging.info(
            'Start vertical interpolation for {} (dimensions={} x {})'.format(myvar, grdROMS.xi_rho, grdROMS.eta_rho))
        outdata = np.empty((outINDEX_ST), dtype=np.float, order='F')

        outdata = interp.interpolation.dovertinter(np.asarray(outdata, order='F'),
                                                   np.asarray(array1, order='F'),
                                                   np.asarray(grdROMS.h, order='F'),
                                                   np.asarray(grdROMS.z_r, order='F'),
                                                   np.asarray(grdMODEL.z_r, order='F'),
                                                   int(grdROMS.nlevels),
                                                   int(grdMODEL.nlevels),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho))

        outdata = np.ma.masked_where(abs(outdata) > 1000, outdata)
        # The BCG has to be capped at 0
        if myvar in ['O3_c', 'O3_TA', 'N1_p', 'N3_p', 'N3_n', 'N5_s', 'O2_o']:
            outdata = np.ma.masked_where(abs(outdata) < 0, outdata)
        # import plotData
        # for k in range(grdROMS.nlevels):
        #     plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, np.squeeze(outdata[k,:,:]),k, varname)

        return outdata

    if myvar == 'vvel':
        logging.info('Start vertical interpolation for uvel (dimensions={} x {})'.format(grdROMS.xi_u, grdROMS.eta_u))
        outdataU = np.zeros((outINDEX_U), dtype=np.float)
        outdataUBAR = np.zeros((outINDEX_UBAR), dtype=np.float)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU, order='F'),
                                                    np.asarray(array1, order='F'),
                                                    np.asarray(grdROMS.h, order='F'),
                                                    np.asarray(grdROMS.z_r, order='F'),
                                                    np.asarray(grdMODEL.z_r, order='F'),
                                                    int(grdROMS.nlevels),
                                                    int(grdMODEL.nlevels),
                                                    int(grdROMS.xi_u),
                                                    int(grdROMS.eta_u),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataU = np.ma.masked_where(abs(outdataU) > 1000, outdataU)

        logging.info('Start vertical interpolation for vvel (dimensions={} x {})'.format(grdROMS.xi_v, grdROMS.eta_v))
        outdataV = np.zeros((outINDEX_V), dtype=np.float)
        outdataVBAR = np.zeros((outINDEX_VBAR), dtype=np.float)

        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV, order='F'),
                                                    np.asarray(array2, order='F'),
                                                    np.asarray(grdROMS.h, order='F'),
                                                    np.asarray(grdROMS.z_r, order='F'),
                                                    np.asarray(grdMODEL.z_r, order='F'),
                                                    int(grdROMS.nlevels),
                                                    int(grdMODEL.nlevels),
                                                    int(grdROMS.xi_v),
                                                    int(grdROMS.eta_v),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataV = np.ma.masked_where(abs(outdataV) > 1000, outdataV)

        z_wu = np.zeros((grdROMS.nlevels + 1, grdROMS.eta_u, grdROMS.xi_u), dtype=np.float)
        z_wv = np.zeros((grdROMS.nlevels + 1, grdROMS.eta_v, grdROMS.xi_v), dtype=np.float)

        outdataUBAR = barotropic.velocity.ubar(np.asarray(outdataU, order='F'),
                                               np.asarray(outdataUBAR, order='F'),
                                               np.asarray(grdROMS.z_w, order='F'),
                                               np.asarray(z_wu, order='F'),
                                               grdROMS.nlevels,
                                               grdROMS.xi_u,
                                               grdROMS.eta_u,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)
        outdataUBAR = np.ma.masked_where(abs(outdataUBAR) > 1000, outdataUBAR)

        # plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataUBAR,1, "ubar")

        outdataVBAR = barotropic.velocity.vbar(np.asarray(outdataV, order='F'),
                                               np.asarray(outdataVBAR, order='F'),
                                               np.asarray(grdROMS.z_w, order='F'),
                                               np.asarray(z_wv, order='F'),
                                               grdROMS.nlevels,
                                               grdROMS.xi_v,
                                               grdROMS.eta_v,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)

        # plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataVBAR,1, "vbar")
        outdataVBAR = np.ma.masked_where(abs(outdataVBAR) > 1000, outdataVBAR)

        return outdataU, outdataV, outdataUBAR, outdataVBAR


def rotate(grdROMS, grdMODEL, data, u, v):
    """
    First rotate the values of U, V at rho points with the angle, and then interpolate
    the rho point values to U and V points and save the result
    """

    urot = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float)
    vrot = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float)

    urot, vrot = interp.interpolation.rotate(np.asarray(urot, order='F'),
                                             np.asarray(vrot, order='F'),
                                             np.asarray(u, order='F'),
                                             np.asarray(v, order='F'),
                                             np.asarray(grdROMS.angle, order='F'),
                                             int(grdROMS.xi_rho),
                                             int(grdROMS.eta_rho),
                                             int(grdMODEL.nlevels))
    return urot, vrot


def interpolate2uv(grdROMS, grdMODEL, urot, vrot):
    Zu = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_u), int(grdROMS.xi_u)), np.float)
    Zv = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_v), int(grdROMS.xi_v)), np.float)

    # Interpolate from RHO points to U and V points for velocities

    Zu = interp.interpolation.rho2u(np.asarray(Zu, order='F'),
                                    np.asarray(urot, order='F'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.nlevels))

    # plotData.contourMap(grdROMS,grdMODEL,Zu[0,:,:],"1",'urot')

    Zv = interp.interpolation.rho2v(np.asarray(Zv, order='F'),
                                    np.asarray(vrot, order='F'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.nlevels))

    # plotData.contourMap(grdROMS,grdMODEL,Zv[0,:,:],"1",'vrot')
    return Zu, Zv


def get_time(confM2R, year, month, day, ntime):
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """

    if confM2R.ocean_indata_type == 'SODA3':
        filename = fc.getSODA3filename(confM2R, year, month, day, None)

    if confM2R.ocean_indata_type == 'SODA3_5DAY':
        filename = fc.getSODA3_5DAYfilename(confM2R, year, month, day, None)

    if confM2R.ocean_indata_type == 'SODAMONTHLY':
        filename = fc.getSODAMONTHLYfilename(confM2R, year, month, None)

    if confM2R.ocean_indata_type == 'GLORYS':
        filename = fc.get_GLORYS_filename(confM2R, year, month, "So")

    if confM2R.ocean_indata_type == 'NORESM':
        filename = fc.getNORESMfilename(confM2R, year, month, "salnlvl")

    # Now open the input file and get the time
    cdf = Dataset(filename)
    jdref = date2num(datetime(1948, 1, 1),
                     units="days since 1948-01-01 00:00:00",
                     calendar="standard")

    if confM2R.ocean_indata_type == 'SODA3_5DAY':
        currentdate = datetime(year, month, day)
        units = confM2R.time_object.units
        jd = date2num(currentdate, units=confM2R.time_object.units, calendar=confM2R.time_object.calendar)

    else:
        # Find the day and month that the GLORYS file represents based on the year and ID number.
        # Each file represents a 1 month average.
        # calendar = cdf.variables["time"].calendar
        units = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, units="days since 1948-01-01 00:00:00", calendar="standard")

    confM2R.grdROMS.time = (jd - jdref)
    confM2R.grdROMS.reftime = jdref
    confM2R.grdROMS.timeunits = "days since 1948-01-01 00:00:00"
    cdf.close()
    logging.info("-------------------------------")
    logging.info('Current time of {} file : {}'.format(confM2R.ocean_indata_type,
                                                         currentdate))
    logging.info("-------------------------------")


def get_3d_data(confM2R, varname, year, month, day, timecounter):
    varN = confM2R.global_varnames.index(varname)

    # The variable splitExtract is defined in IOsubset.py and depends on the orientation
    # and ocean_indata_type of grid (-180-180 or 0-360). Assumes regular grid.

    filename = fc.get_filename(confM2R, year, month, day, confM2R.input_varnames[varN])
    try:
        cdf = Dataset(filename)
    except:
        logging.error("[M2R_model2roms] Unable to open input file {}".format(filename))
        return

    if confM2R.ocean_indata_type == "SODA3":
        data = cdf.variables[confM2R.input_varnames[varN]][month - 1, :, :, :]
        data = np.where(data.mask, confM2R.fillvaluein, data)

    if confM2R.ocean_indata_type == "NORESM":
        # For NorESM data - all data is in one big file so we need the timecounter to access correct data
        myunits = cdf.variables[str(confM2R.input_varnames[varN])].units
        data = np.squeeze(cdf.variables[str(confM2R.input_varnames[varN])][timecounter, :, :, :])
        data = np.where(data.mask, confM2R.fillvaluein, data)

    if confM2R.ocean_indata_type == "GLORYS":
        myunits = cdf.variables[str(confM2R.input_varnames[varN])].units
        data = np.squeeze(cdf.variables[str(confM2R.input_varnames[varN])][0, :, :, :])
        data = np.where(data.mask, confM2R.fillvaluein, data)

    cdf.close()

    if varname == 'temperature' and confM2R.ocean_indata_type in ["GLORYS", "NORESM"]:

        if myunits == "degree_Kelvin" or myunits == "K":
            if confM2R.ocean_indata_type in ["GLORYS"]:
                data = np.where(data <= -32767, confM2R.grdROMS.fillval, data)
            data = data - 273.15

    if confM2R.ocean_indata_type == "GLORYS":
        data = np.where(data <= -32767, confM2R.grdROMS.fillval, data)
        data = np.ma.masked_where(data <= confM2R.grdROMS.fillval, data)

    logging.debug('Data range of {} just after extracting from netcdf file: {:3.3f}-{:3.3f}'.format(
        str(confM2R.input_varnames[varN]),
        float(data.min()), float(data.max())))
    return data


def get_2d_data(confM2R, myvar, year, month, day, timecounter):
    varN = confM2R.global_varnames.index(myvar)

    if confM2R.set_2d_vars_to_zero and confM2R.input_varnames[varN] in ['ageice', 'uice',
                                                                        'vice',
                                                                        'aice',
                                                                        'hice',
                                                                        'hs']:
        return np.zeros((np.shape(confM2R.grdMODEL.lon)))
    else:
        filename = fc.get_filename(confM2R, year, month, day, confM2R.input_varnames[varN])
        try:
            cdf = Dataset(filename)
        except:
            logging.error("[M2R_model2roms] Unable to open input file {}".format(filename))
            return

        if confM2R.ocean_indata_type in ["SODA", "SODA3_5DAY"]:
            data = cdf.variables[confM2R.input_varnames[varN]][0, :, :]

        if confM2R.ocean_indata_type == "SODA3":
            if myvar == 'aice':
                # We only extract the first thickness concentration. Need to fix this so all 5 classes can be extracted.
                # http://www.atmos.umd.edu/~ocean/index_files/soda3_readme.htm
                # hi: sea ice thickness [m ice]
                # mi: sea ice mass [kg/m^2]
                # hs: snow thickness [m snow]
                # {cn1,cn2,cn3,cn4,cn5}: sea ice concentration [0:1] in five ice thickness classes
                data = cdf.variables[confM2R.input_varnames[varN]][int(month - 1), 0, :, :]
            else:
                data = cdf.variables[confM2R.input_varnames[varN]][int(month - 1), :, :]

        if confM2R.ocean_indata_type == "NORESM" and confM2R.set_2d_vars_to_zero is False:
            # myunits = cdf.variables[str(grdROMS.varNames[varN])].units
            # For NORESM data are 12 months of data stored in ice files. Use ID as month indicator to get data.
            data = np.squeeze(cdf.variables[str(confM2R.input_varnames[varN])][timecounter, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if confM2R.ocean_indata_type == "GLORYS":
            data = np.squeeze(cdf.variables[str(confM2R.input_varnames[varN])][0, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if not confM2R.set_2d_vars_to_zero:
            cdf.close()

        if __debug__ and not confM2R.set_2d_vars_to_zero:
            logging.info("[M2R_model2roms] Data range of {} just after extracting from netcdf "
                         "file: {:3.3f}-{:3.3f}".format(str(confM2R.input_varnames[varN]),
                                                        float(data.min()), float(data.max())))
    return data


def convert_MODEL2ROMS(confM2R):
    # First opening of input file is just for initialization of grid
    filenamein = fc.get_filename(confM2R, confM2R.start_year, confM2R.start_month, confM2R.start_day, None)

    # Finalize creating the model grd object now that we know the filename for input data
    confM2R.grdMODEL.create_object(confM2R, filenamein)
    confM2R.grdMODEL.getdims()

    # Create the ESMF weights used to do all of the horizontal interpolation
    interp2D.setup_ESMF_interpolation_weights(confM2R)

    # Now we want to subset the data to avoid storing more information than we need.
    # We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes
    if confM2R.subset_indata:
        IOsubset.find_subset_indices(confM2R.grdMODEL, min_lat=confM2R.subset[0], max_lat=confM2R.subset[1],
                                     min_lon=confM2R.subset[2], max_lon=confM2R.subset[3])

    logging.info("[M2R_model2roms] ==> Initializing done")
    logging.info("[M2R_model2roms] --------------------------")
    logging.info("[M2R_model2roms] ==> Starting loop over time")

    time_counter = 0
    first_run = True

    for year in confM2R.years:
        months = datetimeFunctions.create_list_of_months(confM2R, year)

        for month in months:
            days = datetimeFunctions.create_list_of_days(confM2R, year, month, first_run)

            for day in days:
                # Get the current date for given time-step
                get_time(confM2R, year, month, day, time_counter)

                # Each MODEL file consist only of one time step. Get the subset data selected, and
                # store that time step in a new array:

                if first_run:
                    logging.info("[M2R_model2roms] => NOTE! Make sure that these two arrays are in sequential order:")
                    logging.info("[M2R_model2roms] ==> myvars:     {}".format(confM2R.input_varnames))
                    logging.info("[M2R_model2roms] ==> varNames    {}".format(confM2R.global_varnames))
                    first_run = False

                    if confM2R.subset_indata:
                        # The first iteration we want to organize the subset indices we want to extract
                        # from the input data to get the interpolation correct and to function fast
                        IOsubset.organize_split(confM2R.grdMODEL, confM2R.grdROMS)

                for myvar in confM2R.global_varnames:

                    if myvar in ['temperature', 'salinity', 'uvel', 'vvel', 'O3_c', 'O3_TA', 'N1_p', 'N3_n', 'N5_s',
                                 'O2_o']:
                        data = get_3d_data(confM2R, myvar, year, month, day, time_counter)

                    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
                        data = get_2d_data(confM2R, myvar, year, month, day, time_counter)

                    # Take the input data and horizontally interpolate to your grid
                    array1 = interp2D.do_hor_interpolation_regular_grid(confM2R, data, myvar)

                    if myvar in ['temperature', 'salinity', 'O3_c', 'O3_TA', 'N1_p', 'N3_n', 'N5_s', 'O2_o']:
                        STdata = vertical_interpolation(myvar, array1, array1, confM2R.grdROMS, confM2R.grdMODEL)

                        for dd in range(len(STdata[:, 0, 0])):
                            STdata[dd, :, :] = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval,
                                                        STdata[dd, :, :])

                        STdata = np.where(abs(STdata) > 1000, confM2R.grdROMS.fillval, STdata)

                        IOwrite.write_clim_file(confM2R, time_counter, myvar, STdata)
                        if time_counter == confM2R.grdROMS.inittime and confM2R.grdROMS.write_init is True:
                            IOinitial.create_init_file(confM2R, time_counter, myvar, STdata)

                    if myvar in ['ssh', 'ageice', 'aice', 'hice', 'snow_thick']:
                        SSHdata = array1[0, :, :]

                        SSHdata = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where((abs(SSHdata) > 100) | (SSHdata == 0), confM2R.grdROMS.fillval, SSHdata)

                        # Specific for ROMS - we set 0 where we should have fillvalue for ice otherwise ROMS blows up.
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        IOwrite.write_clim_file(confM2R, time_counter, myvar, SSHdata)

                        if time_counter == confM2R.grdROMS.inittime:
                            IOinitial.create_init_file(confM2R, time_counter, myvar, SSHdata)

                    # The following are special routines used to calculate the u and v velocity
                    # of ice based on the transport, which is divided by snow and ice thickenss
                    # and then multiplied by grid size in dx or dy direction (opposite of transport).
                    if myvar in ['uice', 'vice']:
                        SSHdata = array1[0, :, :]

                        if myvar == "uice":
                            mymask = confM2R.grdROMS.mask_u
                        if myvar == "vice":
                            mymask = confM2R.grdROMS.mask_v

                        SSHdata = np.where(mymask == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where((abs(SSHdata) > 100) | (SSHdata == 0), confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        IOwrite.write_clim_file(confM2R, time_counter, myvar, SSHdata)

                        if time_counter == confM2R.grdROMS.inittime:
                            if myvar in ['uice', 'vice']:
                                IOinitial.create_init_file(confM2R, time_counter, myvar, SSHdata)

                    if myvar == 'uvel':
                        array2 = array1

                    if myvar == 'vvel':
                        urot, vrot = rotate(confM2R.grdROMS, confM2R.grdMODEL, data, array2, array1)
                        u, v = interpolate2uv(confM2R.grdROMS, confM2R.grdMODEL, urot, vrot)

                        Udata, Vdata, UBARdata, VBARdata = vertical_interpolation(myvar, u, v, confM2R.grdROMS,
                                                                                  confM2R.grdMODEL)

                    if myvar == 'vvel':

                        Udata = np.where(confM2R.grdROMS.mask_u == 0, confM2R.grdROMS.fillval, Udata)
                        Udata = np.where(abs(Udata) > 1000, confM2R.grdROMS.fillval, Udata)
                        Vdata = np.where(confM2R.grdROMS.mask_v == 0, confM2R.grdROMS.fillval, Vdata)
                        Vdata = np.where(abs(Vdata) > 1000, confM2R.grdROMS.fillval, Vdata)
                        UBARdata = np.where(confM2R.grdROMS.mask_u == 0, confM2R.grdROMS.fillval, UBARdata)
                        UBARdata = np.where(abs(UBARdata) > 1000, confM2R.grdROMS.fillval, UBARdata)
                        VBARdata = np.where(confM2R.grdROMS.mask_v == 0, confM2R.grdROMS.fillval, VBARdata)
                        VBARdata = np.where(abs(VBARdata) > 1000, confM2R.grdROMS.fillval, VBARdata)

                        IOwrite.write_clim_file(confM2R, time_counter, myvar, Udata, Vdata, UBARdata, VBARdata)

                        if time_counter == confM2R.grdROMS.inittime:
                            IOinitial.create_init_file(confM2R, time_counter, myvar, Udata, Vdata, UBARdata, VBARdata)

                time_counter += 1
