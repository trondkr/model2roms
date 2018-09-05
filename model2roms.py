from __future__ import print_function
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime, timedelta
import numpy as np
import interp2D
import interpolation as interp
import IOwrite
import os
import barotropic
import IOinitial
import IOsubset
import datetimeFunctions

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass
__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2008, 8, 15)
__modified__ = datetime(2014, 12, 1)
__version__ = "1.5"
__status__ = "Development, modified on 15.08.2008,01.10.2009,07.01.2010, 15.07.2014, 01.12.2014, 07.08.2015"


def verticalinterpolation(myvar, array1, array2, grdROMS, grdMODEL):
    outINDEX_ST = (grdROMS.nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    outINDEX_U = (grdROMS.nlevels, grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_V = (grdROMS.nlevels, grdROMS.eta_v, grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v, grdROMS.xi_v)

    if myvar in ['salinity', 'temperature']:
        print('Start vertical interpolation for %s (dimensions=%s x %s)' % (myvar, grdROMS.xi_rho, grdROMS.eta_rho))
        outdata = np.empty((outINDEX_ST), dtype=np.float64, order='Fortran')

        outdata = interp.interpolation.dovertinter(np.asarray(outdata, order='Fortran'),
                                                   np.asarray(array1, order='Fortran'),
                                                   np.asarray(grdROMS.h, order='Fortran'),
                                                   np.asarray(grdROMS.z_r, order='Fortran'),
                                                   np.asarray(grdMODEL.z_r, order='Fortran'),
                                                   int(grdROMS.nlevels),
                                                   int(grdMODEL.nlevels),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho))

        outdata = np.ma.masked_where(abs(outdata) > 1000, outdata)

        # import plotData
        # for k in xrange(len(grdMODEL.h)-1):

        #    plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, np.squeeze(outdata[k,:,:]),k, myvar)

        return outdata

    if myvar == 'vvel':
        print('Start vertical interpolation for uvel (dimensions=%s x %s)' % (grdROMS.xi_u, grdROMS.eta_u))
        outdataU = np.zeros((outINDEX_U), dtype=np.float64)
        outdataUBAR = np.zeros((outINDEX_UBAR), dtype=np.float64)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU, order='Fortran'),
                                                    np.asarray(array1, order='Fortran'),
                                                    np.asarray(grdROMS.h, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.nlevels),
                                                    int(grdMODEL.nlevels),
                                                    int(grdROMS.xi_u),
                                                    int(grdROMS.eta_u),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataU = np.ma.masked_where(abs(outdataU) > 1000, outdataU)

        print('Start vertical interpolation for vvel (dimensions=%s x %s)' % (grdROMS.xi_v, grdROMS.eta_v))
        outdataV = np.zeros((outINDEX_V), dtype=np.float64)
        outdataVBAR = np.zeros((outINDEX_VBAR), dtype=np.float64)

        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV, order='Fortran'),
                                                    np.asarray(array2, order='Fortran'),
                                                    np.asarray(grdROMS.h, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.nlevels),
                                                    int(grdMODEL.nlevels),
                                                    int(grdROMS.xi_v),
                                                    int(grdROMS.eta_v),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataV = np.ma.masked_where(abs(outdataV) > 1000, outdataV)

        z_wu = np.zeros((grdROMS.nlevels + 1, grdROMS.eta_u, grdROMS.xi_u), dtype=np.float64)
        z_wv = np.zeros((grdROMS.nlevels + 1, grdROMS.eta_v, grdROMS.xi_v), dtype=np.float64)

        outdataUBAR = barotropic.velocity.ubar(np.asarray(outdataU, order='Fortran'),
                                               np.asarray(outdataUBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wu, order='Fortran'),
                                               grdROMS.nlevels,
                                               grdROMS.xi_u,
                                               grdROMS.eta_u,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)
        outdataUBAR = np.ma.masked_where(abs(outdataUBAR) > 1000, outdataUBAR)

        # plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataUBAR,1, "ubar")

        outdataVBAR = barotropic.velocity.vbar(np.asarray(outdataV, order='Fortran'),
                                               np.asarray(outdataVBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wv, order='Fortran'),
                                               grdROMS.nlevels,
                                               grdROMS.xi_v,
                                               grdROMS.eta_v,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)

        # plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataVBAR,1, "vbar")
        outdataVBAR = np.ma.masked_where(abs(outdataVBAR) > 1000, outdataVBAR)

        return outdataU, outdataV, outdataUBAR, outdataVBAR


def horizontalinterpolation(confM2R, myvar, data):
    print('Start %s horizontal interpolation for %s' % (confM2R.grdtype, myvar))
    try:
        if myvar in ['temperature', 'salinity']:
            return interp2D.dohorinterpolationregulargrid(confM2R, data)
        elif myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
            return interp2D.dohorinterpolationsshregulargrid(confM2R, data)
        elif myvar in ['uvel', 'vvel']:
            return interp2D.dohorinterpolationregulargrid(confM2R, data)
    except IOError as error:
        print("An error occurred in horizontalinterpolation: {}".format(error))
        raise


def rotate(grdROMS, grdMODEL, data, u, v):
    """
    First rotate the values of U, V at rho points with the angle, and then interpolate
    the rho point values to U and V points and save the result
    """

    urot = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)
    vrot = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)

    urot, vrot = interp.interpolation.rotate(np.asarray(urot, order='Fortran'),
                                             np.asarray(vrot, order='Fortran'),
                                             np.asarray(u, order='Fortran'),
                                             np.asarray(v, order='Fortran'),
                                             np.asarray(grdROMS.angle, order='Fortran'),
                                             int(grdROMS.xi_rho),
                                             int(grdROMS.eta_rho),
                                             int(grdMODEL.nlevels))
    return urot, vrot


def interpolate2uv(grdROMS, grdMODEL, urot, vrot):
    Zu = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_u), int(grdROMS.xi_u)), np.float64)
    Zv = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_v), int(grdROMS.xi_v)), np.float64)

    # Interpolate from RHO points to U and V points for velocities

    Zu = interp.interpolation.rho2u(np.asarray(Zu, order='Fortran'),
                                    np.asarray(urot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.nlevels))

    # plotData.contourMap(grdROMS,grdMODEL,Zu[0,:,:],"1",'urot')

    Zv = interp.interpolation.rho2v(np.asarray(Zv, order='Fortran'),
                                    np.asarray(vrot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.nlevels))

    # plotData.contourMap(grdROMS,grdMODEL,Zv[0,:,:],"1",'vrot')

    return Zu, Zv


def getTime(confM2R, year, month, day):
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    if confM2R.indatatype == 'SODA':
        filename = getSODAfilename(confM2R, year, month, None)

    if confM2R.indatatype == 'SODA3':
        filename = getSODA3filename(confM2R, year, month, None)

    if confM2R.indatatype == 'SODAMONTHLY':
        filename = getSODAMONTHLYfilename(confM2R, year, month, None)

    if confM2R.indatatype == 'GLORYS':
        filename = getGLORYSfilename(confM2R, year, month, "S")

    if confM2R.indatatype == 'WOAMONTHLY':
        filename = getWOAMONTHLYfilename(confM2R, year, month, "temperature")

    if confM2R.indatatype == 'NORESM':
        filename = getNORESMfilename(confM2R, year, month, "saln")

    if confM2R.indatatype == 'NS8KM':
        filename = getNS8KMfilename(confM2R, year, month, "salt")

    if confM2R.indatatype == 'NS8KMZ':
        filename, readFromOneFile = getNS8KMZfilename(confM2R, year, month, "salt")

    # Now open the input file and get the time
    cdf = Dataset(filename)

    if confM2R.indatatype == 'NORESM':
        jdref = date2num(datetime(1800, 1, 1), cdf.variables["time"].units, calendar=cdf.variables["time"].calendar)
    elif confM2R.indatatype == 'NS8KMZ':
        jdref = date2num(datetime(1948, 1, 1), units="days since 1948-01-01 00:00:00", calendar="standard")
    elif confM2R.indatatype == 'GLORYS':
        jdref = date2num(datetime(1948, 1, 1), cdf.variables["time_counter"].units,
                         calendar=cdf.variables["time_counter"].calendar)
    elif confM2R.indatatype == 'NS8KM':
        jdref = date2num(datetime(1948, 1, 1), cdf.variables["ocean_time"].units,
                         calendar=cdf.variables["ocean_time"].calendar)
    elif confM2R.indatatype == 'SODA3':
        jdref = date2num(datetime(1948, 1, 1), units="days since 1948-01-01 00:00:00", calendar="standard")
    else:
        jdref = date2num(datetime(1948, 1, 1), cdf.variables["time"].units, calendar=cdf.variables["time"].calendar)

    if confM2R.indatatype == 'SODAMONTHLY':
        # Find the day and month that the SODAMONTHLY file respresents based on the year and ID number.
        # Each SODA file represents a 1 month average.

        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if confM2R.indatatype == 'SODA3':
        # Each SODA file represents 12 month averages.

        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, units="days since 1948-01-01 00:00:00", calendar="standard")

    if confM2R.indatatype == 'GLORYS':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
        # Each file represents a 1 month average.
        mycalendar = cdf.variables["time_counter"].calendar
        myunits = cdf.variables["time_counter"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if confM2R.indatatype == 'NS8KM':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
        # Each file represents a 1 month average.
        mycalendar = cdf.variables["ocean_time"].calendar
        myunits = cdf.variables["ocean_time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if confM2R.indatatype == 'NS8KMZ':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
        # Each file represents a 1 month average.
        currentdate = datetime(year, month, day)
        myunits = cdf.variables["time"].units
        jd = date2num(currentdate, myunits, calendar="gregorian")
        print("Days:", jd, currentdate, year, month, day)

    if confM2R.indatatype == 'NORESM':
        # Find the day and month that the NORESM file. We need to use the time modules from
        # netcdf4 for python as they handle calendars that are no_leap.
        # http://www.esrl.noaa.gov/psd/people/jeffrey.s.whitaker/python/netcdftime.html#datetime
        mydays = cdf.variables["time"][0]
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = num2date(mydays, units=myunits, calendar=mycalendar)
        jd = date2num(currentdate, myunits, calendar='noleap')

    confM2R.grdROMS.time = (jd - jdref)
    confM2R.grdROMS.reftime = jdref
    confM2R.grdROMS.timeunits = myunits
    cdf.close()
    print("-------------------------------")
    print('\nCurrent time of %s file : %s' % (confM2R.indatatype, currentdate))
    print("-------------------------------")

def getGLORYSfilename(confM2R, year, month, myvar):
    # Month indicates month
    # myvar:S,T,U,V
    print("MYVAR: %s - %s" % (myvar, confM2R.modelpath))
    if myvar in ['iicevelu', 'iicevelv', 'ileadfra', 'iicethic']:
        myvarPrefix = 'icemod'
        myvar = "ice"
    elif myvar in ['sossheig']:

        if (confM2R.start_year == 2014):
            myvarPrefix = 'grid2D'
        elif (2013 > confM2R.start_year >= 2010 and confM2R.start_month == 12) or (2013 > confM2R.start_year >= 2011):
            myvarPrefix = 'SSH'
        elif confM2R.start_year == 2013:
            myvarPrefix = 'SSH'
        myvar = "ssh"
    elif (myvar in ['vozocrtx', 'vomecrty']):
        myvar = 'u-v'
        myvarPrefix = 'gridUV'
    elif myvar in ['votemper']:
        myvar = 't'
        myvarPrefix = 'gridT'
    elif myvar in ['vosaline']:
        myvar = 's'
        myvarPrefix = 'gridS'
    else:
        myvarPrefix = 'grid' + str(myvar.upper())

    # GLORYS change the name in the middle of the time-series (on December 2010) and we have to
    # account for that
    if confM2R.start_year == 2014:
        production = "R20151218"
    elif (confM2R.start_2013 > confM2R.start_year >= 2010 and confM2R.start_month == 12) or (
            2013 > confM2R.start_year >= 2011):
        production = "R20140520"
    elif confM2R.start_year == 2013:
        production = "R20141205"
    else:
        production = "R20130808"

    if confM2R.start_month < 10:
        filename = confM2R.modelpath + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' \
                   + str(myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(confM2R.start_year) + '0' + str(
            confM2R.start_month) \
                   + '15_' + str(production) + '_' + str(myvarPrefix) + '.nc'

    if confM2R.start_month >= 10:
        filename = confM2R.modelpath + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' \
                   + str(myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(confM2R.start_year) + str(confM2R.start_month) \
                   + '15_' + str(production) + '_' + str(myvarPrefix) + '.nc'

    return filename


def getNORESMfilename(confM2R, year, month, myvar):
    if myvar == 'grid':
        filename = confM2R.modelpath + 'GRID/NorESM.nc'
    else:
        if myvar in ['iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']:
            filename = confM2R.modelpath + 'ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h.' + str(year) + '.nc'
        else:
            if confM2R.start_month < 10:
                filename = confM2R.modelpath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.' + str(year) + '-0' + str(
                    ID) + '.nc'
            else:
                filename = confM2R.modelpath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.' + str(year) + '-' + str(
                    ID) + '.nc'

    return filename


def getNS8KMfilename(confM2R, year, month, myvar):
    if (confM2R.start_month < 10):
        mymonth = '0%s' % (confM2R.start_month)
    else:
        mymonth = '%s' % (confM2R.start_month)
    filename = confM2R.modelpath + str(confM2R.start_year) + str(
        mymonth) + '15_mm-IMR-MODEL-ROMS-NWS-20140430-fv02.1.nc'

    return filename, mymonth


def getNS8KMZfilename(confM2R, year, month, myvar):
    # allInOneFile = '/work/users/trondk/KINO/FORCING/1600M/northsea_8km_z_mean.nc_2010-2013.nc'
    allInOneFile = ''

    if os.path.exists(allInOneFile):
        readFromOneFile = True
        print("NOTE ! READING ALL MYOCEAN FORCING DATA FROM ONE FILE")
        return allInOneFile, readFromOneFile
    else:
        readFromOneFile = False
        if (confM2R.start_month < 10):
            mymonth = '0%s' % (confM2R.start_month)
        else:
            mymonth = '%s' % (confM2R.start_month)
        filename = confM2R.modelpath + str(confM2R.start_year) + str(
            mymonth) + '15_mm-IMR-MODEL-ROMS-NWS-20160203-fv02.1.nc'

        return filename, readFromOneFile


def getSODAMONTHLYfilename(confM2R, year, month, myvar):
    if confM2R.start_month < 10:
        return confM2R.modelpath + 'SODA_2.0.2_' + str(year) + '0' + str(month) + '.cdf'
    elif confM2R.start_month >= 10:
        return confM2R.modelpath + 'SODA_2.0.2_' + str(year) + str(month) + '.cdf'


def getSODAfilename(confM2R, year, month, myvar):
    return confM2R.modelpath + "SODA_2.0.2_" + str(year) + "_" + str(month) + ".cdf"


def getSODA3filename(confM2R, year, month, myvar):
    if (myvar in ['cn', 'hi', 'hs']):
        return confM2R.modelpath + "soda3.3.1_mn_ice_reg_" + str(year) + ".nc"
    else:
        return confM2R.modelpath + "soda3.3.1_mn_ocean_reg_" + str(year) + ".nc"


def getWOAMONTHLYfilename(confM2R, year, month, myvar):
    if myvar == "temperature":
        return confM2R.modelpath + 'temperature_monthly_1deg.nc'
    elif myvar == "salinity":
        return confM2R.modelpath + 'salinity_monthly_1deg.nc'
    else:
        print("Could not find any input files in folder: %s" % confM2R.modelpath)


def get3ddata(confM2R, myvar, year, month, day):
    if myvar == 'temperature':
        varN = 0
    if myvar == 'salinity':
        varN = 1
    if myvar == 'uvel':
        varN = 3
    if myvar == 'vvel':
        varN = 4

    # The variable splitExtract is defined in IOsubset.py and depends on the orientation
    # and indatatype of grid (-180-180 or 0-360). Assumes regular grid.
    if confM2R.useesmf:
        if confM2R.indatatype == "SODA":
            filename = getSODAfilename(confM2R, year, month, None)
            cdf = Dataset(filename)
            data = cdf.variables[confM2R.inputdatavarnames[varN]][0, :, :, :]

        if confM2R.indatatype == "SODA3":
            filename = getSODA3filename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            cdf = Dataset(filename)
            print("=>Extracting data for month %s from SODA3 %s " % (month - 1, filename))
            data = cdf.variables[confM2R.inputdatavarnames[varN]][month - 1, :, :, :]

        if confM2R.indatatype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            cdf = Dataset(filename)
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][:, :, :]

        if confM2R.indatatype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            cdf = Dataset(filename)
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][month - 1, :, :, :]

        if confM2R.indatatype == "NORESM":
            cdf = Dataset(getNORESMfilename(confM2R, year, month, confM2R.inputdatavarnames[varN]))
            myunits = cdf.variables[str(confM2R.inputdatavarnames[varN])].units
            data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

            print("Data range", np.min(data), np.max(data))

        if confM2R.indatatype == "NS8KMZ":
            filename, readFromOneFile = getNS8KMZfilename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            cdf = Dataset(filename)
            print("Reading from one file %s" % (readFromOneFile))

            myunits = cdf.variables[str(confM2R.inputdatavarnames[varN])].units
            if (readFromOneFile):
                currentdate = datetime(year, month, day, 12)
                jd = date2num(currentdate, cdf.variables["time"].units, calendar="gregorian")
                print("Days:", jd, currentdate)

                timesteps = (cdf.variables["time"][:]).tolist()
                timeindex = timesteps.index(jd)
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][timeindex, :, :, :])
            else:
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :, :])
                print("Range of data", varN, np.min(data), np.max(data))
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if confM2R.indatatype == "GLORYS":
            cdf = Dataset(getGLORYSfilename(confM2R, year, month, confM2R.inputdatavarnames[varN]))
            myunits = cdf.variables[str(confM2R.inputdatavarnames[varN])].units
            data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        cdf.close()

    if myvar == 'temperature' and confM2R.indatatype in ["NS8KMZ", "GLORYS", "NORESM"]:

        if myunits == "degree_Kelvin" or myunits == "K":
            if confM2R.indatatype in ["GLORYS"]:
                data = np.where(data <= -32.767, confM2R.grdROMS.fillval, data)
            data = data - 273.15

    if confM2R.indatatype == "GLORYS":
        data = np.where(data <= -32.767, confM2R.grdROMS.fillval, data)
        data = np.ma.masked_where(data <= confM2R.grdROMS.fillval, data)

    if __debug__:
        print("Data range of %s just after extracting from netcdf file: %s - %s" % (str(confM2R.inputdatavarnames[varN]),
                                                                                    data.min(), data.max()))

    return data


def get2ddata(confM2R, myvar, year, month, day):
    indexROMS_SSH = (confM2R.grdROMS.eta_rho, confM2R.grdROMS.xi_rho)

    if confM2R.indatatype == "NORESM":
        if myvar == 'ageice':
            varN = 5
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'uice':
            varN = 6
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'vice':
            varN = 7
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'aice':
            varN = 8
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':
            varN = 9
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'snow_thick':
            varN = 10
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if confM2R.indatatype == "GLORYS":
        if myvar == 'uice':
            varN = 5
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'vice':
            varN = 6
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'aice':
            varN = 7
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':
            varN = 8
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    print("Current type %s and variable %s" % (confM2R.indatatype, myvar))
    if confM2R.indatatype == "SODA3":
        if myvar == 'aice':
            varN = 5
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':
            varN = 6
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'snow_thick':
            varN = 7
            SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if myvar == 'ssh':
        varN = 2
        SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if confM2R.useesmf:
        if confM2R.indatatype == "SODA":
            filename = getSODAfilename(confM2R, year, month, day, None)
            cdf = Dataset(filename)
            data = cdf.variables[confM2R.inputdatavarnames[varN]][0, :, :]

        if confM2R.indatatype == "SODA3":
            filename = getSODA3filename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            print("Trying to open file %s" % (filename))
            cdf = Dataset(filename)
            if myvar == 'aice':
                # We only extract the first thickness concentration. Need to fix this so all 5 classes can be extracted.
                # http://www.atmos.umd.edu/~ocean/index_files/soda3_readme.htm
                # hi: sea ice thickness [m ice]
                # mi: sea ice mass [kg/m^2]
                # hs: snow thickness [m snow]
                # {cn1,cn2,cn3,cn4,cn5}: sea ice concentration [0:1] in five ice thickness classes
                data = cdf.variables[confM2R.inputdatavarnames[varN]][int(month - 1), 0, :, :]
            else:
                data = cdf.variables[confM2R.inputdatavarnames[varN]][int(month - 1), :, :]

        if confM2R.indatatype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(confM2R, year, month, None)
            cdf = Dataset(filename)
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][:, :]

        if confM2R.indatatype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(confM2R, year, month, myvar)
            cdf = Dataset(filename)
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][month - 1, :, :]

        if confM2R.indatatype == "NORESM":
            cdf = Dataset(getNORESMfilename(confM2R, year, month, confM2R.inputdatavarnames[varN]))
            # myunits = cdf.variables[str(grdROMS.varNames[varN])].units
            # For NORESM data are 12 months of data stored in ice files. Use ID as month indicator to get data.
            if myvar in ['ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][month - 1, :, :])
            else:
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if confM2R.indatatype == "GLORYS":
            cdf = Dataset(getGLORYSfilename(confM2R, year, month, confM2R.inputdatavarnames[varN]))
            data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if confM2R.indatatype == "NS8KMZ":
            filename, readFromOneFile = getNS8KMZfilename(confM2R, year, month, confM2R.inputdatavarnames[varN])
            cdf = Dataset(filename)
        
            if (readFromOneFile):
                jdref = date2num(datetime(1948, 1, 1), cdf.variables["time"].units,
                                 calendar=cdf.variables["time"].calendar)
                currentdate = datetime(year, month, day, 12)
                jd = date2num(currentdate, cdf.variables["time"].units, calendar="gregorian")

                timesteps = (cdf.variables["time"][:]).tolist()
                timeindex = timesteps.index(jd)
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][timeindex, :, :])
            else:
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0, :, :])

            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

            print("Extracted raw data: %s min: %s max: %s" % (myvar, np.min(data), np.max(data)))
        cdf.close()

        if __debug__:
            print("Data range of %s just after extracting from netcdf file: %s - %s" % (str(grdROMS.varNames[varN]),
                                                                                        data.min(), data.max()))

    return data


def convertMODEL2ROMS(confM2R):
    # First opening of input file is just for initialization of grid
    if confM2R.indatatype == 'SODA':
        filenamein = getSODAfilename(confM2R, confM2R.start_year, confM2R.start_month, "salinity")
    if confM2R.indatatype == 'SODA3':
        filenamein = getSODA3filename(confM2R, confM2R.start_year, confM2R.start_month, "salinity")
    if confM2R.indatatype == 'SODAMONTHLY':
        filenamein = getSODAfilename(confM2R, confM2R.start_year, confM2R.start_month, "salinity")
    if confM2R.indatatype == 'NORESM':
        filenamein = getNORESMfilename(confM2R, confM2R.start_year, confM2R.start_month, "grid")
    if confM2R.indatatype == 'WOAMONTHLY':
        filenamein = getWOAMONTHLYfilename(confM2R, confM2R.start_year, confM2R.start_month, "temperature")
    if confM2R.indatatype == 'GLORYS':
        filenamein = getGLORYSfilename(confM2R, confM2R.start_year, confM2R.start_month, "S")
    if confM2R.indatatype == 'GLORYS':
        filenamein = getGLORYSfilename(confM2R, confM2R.start_year, confM2R.start_month, "S")
    if confM2R.indatatype == 'NS8KM':
        filenamein = getNS8KMfilename(confM2R, confM2R.start_year, confM2R.start_month, "S")
    if confM2R.indatatype == 'NS8KMZ':
        filenamein, readFromOneFile = getNS8KMZfilename(confM2R, confM2R.start_year, confM2R.start_month, "S")

    # Finalize creating the model grd object now that we know the filename for input data
    confM2R.grdMODEL.opennetcdf(filenamein)
    confM2R.grdMODEL.createobject(confM2R)
    confM2R.grdMODEL.getdims()

    if confM2R.useesmf:
        print("=>Creating the interpolation weights and indexes using ESMF (this may take some time....):")

        print("  -> regridSrc2Dst at RHO points")
        confM2R.grdMODEL.fieldSrc = ESMF.Field(confM2R.grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_rho = ESMF.Field(confM2R.grdROMS.esmfgrid, "fieldDst",
                                                   staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_rho,
                                                         regrid_method=ESMF.RegridMethod.BILINEAR,
                                                         unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("  -> regridSrc2Dst at U points")
        confM2R.grdMODEL.fieldSrc = ESMF.Field(confM2R.grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_u = ESMF.Field(confM2R.grdROMS.esmfgrid_u, "fieldDst_u",
                                                 staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_u = ESMF.Regrid(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_u,
                                                       regrid_method=ESMF.RegridMethod.BILINEAR,
                                                       unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("  -> regridSrc2Dst at V points")
        confM2R.grdMODEL.fieldSrc = ESMF.Field(confM2R.grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_v = ESMF.Field(confM2R.grdROMS.esmfgrid_v, "fieldDst_v",
                                                 staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_v = ESMF.Regrid(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_v,
                                                       regrid_method=ESMF.RegridMethod.BILINEAR,
                                                       unmapped_action=ESMF.UnmappedAction.IGNORE)

    # Now we want to subset the data to avoid storing more information than we need.
    # We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes
    if confM2R.subsetindata:
        IOsubset.findSubsetIndices(confM2R.grdMODEL, min_lat=confM2R.subset[0], max_lat=confM2R.subset[1],
                                   min_lon=confM2R.subset[2], max_lon=confM2R.subset[3])

    print('==> Initializing done')
    print('\n--------------------------')
    print('==> Starting loop over time')

    time = 0
    firstrun = True

    for year in confM2R.years:
        months = datetimeFunctions.createlistofmonths(confM2R, year)

        for month in months:
            days = datetimeFunctions.createlistofdays(confM2R, year, month)

            for day in days:
                # Get the current date for given timestep 
                getTime(confM2R, year, month, day)

                # Each MODEL file consist only of one time step. Get the subset data selected, and
                # store that time step in a new array:

                if firstrun:
                    print("=> NOTE! Make sure that these two arrays are in sequential order:")
                    print("==> myvars:     %s" % confM2R.inputdatavarnames)
                    print("==> varNames    %s" % confM2R.globalvarnames)
                    firstrun = False

                    if confM2R.subsetindata:
                        # The first iteration we want to organize the subset indices we want to extract
                        # from the input data to get the interpolation correct and to function fast
                        IOsubset.organizeSplit(confM2R.grdMODEL, confM2R.grdROMS)

                for myvar in confM2R.globalvarnames:

                    if myvar in ['temperature', 'salinity', 'uvel', 'vvel']:
                        data = get3ddata(confM2R, myvar, year, month, day)

                    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
                        data = get2ddata(confM2R, myvar, year, month, day)

                    # Take the input data and horizontally interpolate to your grid

                    array1 = horizontalinterpolation(confM2R, myvar, data)

                    if myvar in ['temperature', 'salinity']:
                        STdata = verticalinterpolation(myvar, array1, array1, confM2R.grdROMS, confM2R.grdMODEL)

                        for dd in range(len(STdata[:, 0, 0])):
                            STdata[dd, :, :] = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval,
                                                        STdata[dd, :, :])

                        STdata = np.where(abs(STdata) > 1000, confM2R.grdROMS.fillval, STdata)

                        IOwrite.writeclimfile(confM2R, time, myvar, STdata)
                        if time == confM2R.grdROMS.inittime and confM2R.grdROMS.write_init is True:
                            IOinitial.createinitfile(confM2R, time, myvar, STdata)

                    if myvar in ['ssh', 'ageice', 'aice', 'hice', 'snow_thick']:
                        SSHdata = array1[0, :, :]

                        SSHdata = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, confM2R.grdROMS.fillval, SSHdata)

                        # Specific for ROMs. We set 0 where we should have fillvalue for ice otherwise ROMS blows up.
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        IOwrite.writeclimfile(confM2R, time,  myvar, SSHdata)

                        if time == confM2R.grdROMS.inittime:
                            IOinitial.createinitfile(confM2R, time,  myvar, SSHdata)

                    # The following are special routines used to calculate the u and v velocity
                    # of ice based on the transport, which is divided by snow and ice thickenss
                    # and then multiplied by grid size in dx or dy direction (opposite of transport).
                    if myvar in ['uice', 'vice']:
                        SSHdata = array1[0, :, :]

                        if myvar == "uice": mymask = confM2R.grdROMS.mask_u
                        if myvar == "vice": mymask = confM2R.grdROMS.mask_v

                        SSHdata = np.where(mymask == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        # SSHdata = np.ma.masked_where(abs(SSHdata) > 1000, SSHdata)

                        print("Data range of %s after interpolation: %3.3f to %3.3f" % (
                        myvar, SSHdata.min(), SSHdata.max()))

                        IOwrite.writeclimfile(confM2R, time, myvar, SSHdata)

                        if time == confM2R.grdROMS.inittime:
                            if myvar == 'uice':
                                IOinitial.createinitfile(confM2R, time,  myvar, SSHdata)
                            if myvar == 'vice':
                                IOinitial.createinitfile(confM2R, time,  myvar, SSHdata)

                    if myvar == 'uvel':
                        array2 = array1

                    if myvar == 'vvel':
                        indexROMS_UBAR = (confM2R.grdROMS.eta_u, confM2R.grdROMS.xi_u)
                        indexROMS_VBAR = (confM2R.grdROMS.eta_v, confM2R.grdROMS.xi_v)
                      #  UBARdata = np.zeros((indexROMS_UBAR), dtype=np.float64)
                      #  VBARdata = np.zeros((indexROMS_VBAR), dtype=np.float64)

                        urot, vrot = rotate(confM2R.grdROMS, confM2R.grdMODEL, data, array2, array1)

                        u, v = interpolate2uv(confM2R.grdROMS, confM2R.grdMODEL, urot, vrot)

                        Udata, Vdata, UBARdata, VBARdata = verticalinterpolation(myvar, u, v, confM2R.grdROMS,
                                                                                 confM2R.grdMODEL)

                    if myvar == 'vvel':
                        #   print "Data range of U after interpolation: %3.3f to %3.3f - V after scaling: %3.3f to %3.3f" % (
                        #      Udata.min(), Udata.max(), Vdata.min(), Vdata.max())

                        Udata = np.where(confM2R.grdROMS.mask_u == 0, confM2R.grdROMS.fillval, Udata)
                        Udata = np.where(abs(Udata) > 1000, confM2R.grdROMS.fillval, Udata)
                        Vdata = np.where(confM2R.grdROMS.mask_v == 0, confM2R.grdROMS.fillval, Vdata)
                        Vdata = np.where(abs(Vdata) > 1000, confM2R.grdROMS.fillval, Vdata)
                        UBARdata = np.where(confM2R.grdROMS.mask_u == 0, confM2R.grdROMS.fillval, UBARdata)
                        UBARdata = np.where(abs(UBARdata) > 1000, confM2R.grdROMS.fillval, UBARdata)
                        VBARdata = np.where(confM2R.grdROMS.mask_v == 0, confM2R.grdROMS.fillval, VBARdata)
                        VBARdata = np.where(abs(VBARdata) > 1000, confM2R.grdROMS.fillval, VBARdata)

                        IOwrite.writeclimfile(confM2R, time, myvar,  Udata, Vdata, UBARdata, VBARdata)

                        if time == confM2R.grdROMS.inittime:
                            IOinitial.createinitfile(confM2R, time, myvar, Udata, Vdata, UBARdata, VBARdata)

                time += 1
