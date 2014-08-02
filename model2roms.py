from netCDF4 import Dataset, datetime, date2num,num2date
import numpy as np
import interp2D
import interpolation as interp
import IOwrite

import grd
import barotropic
import IOinitial
import IOsubset

try:
    import ESMF
except ImportError:
    print "Could not find module ESMF"
    pass
__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@imr.no'
__created__ = datetime(2008, 8, 15)
__modified__ = datetime(2014, 7, 15)
__version__ = "1.5"
__status__ = "Development, modified on 15.08.2008,01.10.2009,07.01.2010, 15.07.2014"


def VerticalInterpolation(myvar, array1, array2, grdROMS, grdMODEL):
    outINDEX_ST = (grdROMS.Nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    outINDEX_U = (grdROMS.Nlevels, grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_V = (grdROMS.Nlevels, grdROMS.eta_v, grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v, grdROMS.xi_v)

    if myvar in ['salinity', 'temperature']:
        print 'Start vertical interpolation for %s (dimensions=%s x %s)' % (myvar, grdROMS.xi_rho, grdROMS.eta_rho)
        outdata = np.zeros((outINDEX_ST), dtype=np.float64, order='Fortran')

        outdata = interp.interpolation.dovertinter(np.asarray(outdata, order='Fortran'),
                                                   np.asarray(array1, order='Fortran'),
                                                   np.asarray(grdROMS.depth, order='Fortran'),
                                                   np.asarray(grdROMS.z_r, order='Fortran'),
                                                   np.asarray(grdMODEL.z_r, order='Fortran'),
                                                   int(grdROMS.Nlevels),
                                                   int(grdMODEL.Nlevels),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho))

        outdata = np.ma.masked_where(abs(outdata) > 1000, outdata)


        #import plotData
        #for k in range(grdROMS.Nlevels):
        #plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, np.squeeze(outdata[34,:,:]),34, myvar)


        return outdata

    if myvar == 'vvel':
        print 'Start vertical interpolation for uvel (dimensions=%s x %s)' % (grdROMS.xi_u, grdROMS.eta_u)
        outdataU = np.zeros((outINDEX_U), dtype=np.float64)
        outdataUBAR = np.zeros((outINDEX_UBAR), dtype=np.float64)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU, order='Fortran'),
                                                    np.asarray(array1, order='Fortran'),
                                                    np.asarray(grdROMS.depth, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.Nlevels),
                                                    int(grdMODEL.Nlevels),
                                                    int(grdROMS.xi_u),
                                                    int(grdROMS.eta_u),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataU = np.ma.masked_where(abs(outdataU) > 1000, outdataU)

        print 'Start vertical interpolation for vvel (dimensions=%s x %s)' % (grdROMS.xi_v, grdROMS.eta_v)
        outdataV = np.zeros((outINDEX_V), dtype=np.float64)
        outdataVBAR = np.zeros((outINDEX_VBAR), dtype=np.float64)

        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV, order='Fortran'),
                                                    np.asarray(array2, order='Fortran'),
                                                    np.asarray(grdROMS.depth, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.Nlevels),
                                                    int(grdMODEL.Nlevels),
                                                    int(grdROMS.xi_v),
                                                    int(grdROMS.eta_v),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataV = np.ma.masked_where(abs(outdataV) > 1000, outdataV)

        z_wu = np.zeros((grdROMS.Nlevels + 1, grdROMS.eta_u, grdROMS.xi_u), dtype=np.float64)
        z_wv = np.zeros((grdROMS.Nlevels + 1, grdROMS.eta_v, grdROMS.xi_v), dtype=np.float64)

        outdataUBAR = barotropic.velocity.ubar(np.asarray(outdataU, order='Fortran'),
                                               np.asarray(outdataUBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wu, order='Fortran'),
                                               grdROMS.Nlevels,
                                               grdROMS.xi_u,
                                               grdROMS.eta_u,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)
        outdataUBAR = np.ma.masked_where(abs(outdataUBAR) > 1000, outdataUBAR)

        #plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataUBAR,1, "ubar")



        outdataVBAR = barotropic.velocity.vbar(np.asarray(outdataV, order='Fortran'),
                                               np.asarray(outdataVBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wv, order='Fortran'),
                                               grdROMS.Nlevels,
                                               grdROMS.xi_v,
                                               grdROMS.eta_v,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)

        #plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataVBAR,1, "vbar")
        outdataVBAR = np.ma.masked_where(abs(outdataVBAR) > 1000, outdataVBAR)

        return outdataU, outdataV, outdataUBAR, outdataVBAR


def HorizontalInterpolation(myvar, grdROMS, grdMODEL, data, show_progress):
    print 'Start %s horizontal interpolation for %s' % (grdMODEL.grdType, myvar)

    if myvar in ['temperature', 'salinity']:
        array1 = interp2D.doHorInterpolationRegularGrid(myvar, grdROMS, grdMODEL, data, show_progress)
    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
        array1 = interp2D.doHorInterpolationSSHRegularGrid(myvar, grdROMS, grdMODEL, data)
    if myvar in ['uvel', 'vvel']:
        array1 = interp2D.doHorInterpolationRegularGrid(myvar, grdROMS, grdMODEL, data, show_progress)

    return array1


def rotate(grdROMS, grdMODEL, data, u, v):
    """
    First rotate the values of U, V at rho points with the angle, and then interpolate
    the rho point values to U and V points and save the result
    """

    urot = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)
    vrot = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)

    urot, vrot = interp.interpolation.rotate(np.asarray(urot, order='Fortran'),
                                             np.asarray(vrot, order='Fortran'),
                                             np.asarray(u, order='Fortran'),
                                             np.asarray(v, order='Fortran'),
                                             np.asarray(grdROMS.angle, order='Fortran'),
                                             int(grdROMS.xi_rho),
                                             int(grdROMS.eta_rho),
                                             int(grdMODEL.Nlevels))
    return urot, vrot


def interpolate2UV(grdROMS, grdMODEL, urot, vrot):
    Zu = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_u), int(grdROMS.xi_u)), np.float64)
    Zv = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_v), int(grdROMS.xi_v)), np.float64)

    # Interpolate from RHO points to U and V points for velocities

    Zu = interp.interpolation.rho2u(np.asarray(Zu, order='Fortran'),
                                    np.asarray(urot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.Nlevels))


    #plotData.contourMap(grdROMS,grdMODEL,Zu[0,:,:],"1",'urot')

    Zv = interp.interpolation.rho2v(np.asarray(Zv, order='Fortran'),
                                    np.asarray(vrot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.Nlevels))


    #plotData.contourMap(grdROMS,grdMODEL,Zv[0,:,:],"1",'vrot')

    return Zu, Zv


def getTime(cdf, grdROMS, grdMODEL, year, ID, mytime, mytype, firstRun):
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """

    if (mytype)=='NORESM':
       jdref = date2num(datetime(1800,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)
    else:
        jdref = date2num(datetime(1948,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)

    if mytype == 'SODA':

        # Find the day and month that the SODA file respresents based on the year and ID number.
        # Each SODA file represents a 5 day average, therefore we let the date we find be the first day
        # of those 5 days. Thats the reason we subtract 4 below for day of month.
        import date

        days = 0.0; month = 1; loop = True

        while loop is True:
            d = date.NumberDaysMonth(month, year)
            if days + d < int(ID) * 5:
                days = days + d
                month += 1
            else:
                day = int(int(ID) * 5 - days)
                loop = False

        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if mytype == 'SODAMONTHLY':
        # Find the day and month that the SODAMONTHLY file respresents based on the year and ID number.
        # Each SODA file represents a 1 month average.

        month = ID
        day = 15
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)


    if mytype == 'GLORYS2V1':
        # Find the day and month that the SODAMONTHLY file respresents based on the year and ID number.
        # Each SODA file represents a 1 month average.

        month = ID
        day = 15
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)


    if mytype == 'NORESM':
        # Find the day and month that the NORESM file. We need to use the time modules from
        # netcdf4 for python as they handle calendars that are no_leap.
        # http://www.esrl.noaa.gov/psd/people/jeffrey.s.whitaker/python/netcdftime.html#datetime
        mydays = cdf.variables["time"][0]
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        # Fake the start date first time around
        if (firstRun):
            currentdate = datetime(2006,1,1)
            print "NOTICE!\n First datestamp in result files are hardcoded to %s"%(currentdate)
        else:
            currentdate = num2date(mydays, units=myunits, calendar=mycalendar)
        jd = date2num(currentdate, myunits, calendar='noleap')

    grdROMS.time = (jd - jdref)
    grdROMS.reftime = jdref

    print '\nCurrent time of %s file : %s' % (mytype, currentdate)


def getGLORYS2V1filename(year, ID, myvar, dataPath):
    # Month indicates month
    # myvar:S,T,U,V
    if ID < 10: filename = dataPath + 'global-reanalysis-phys-001-009-ran-fr-glorys2-grid' + str(
        myvar.lower()) + '/REGRID/GLORYS2V1_ORCA025_' + str(year) + '0' + str(ID) + '15_R20110216_grid' + str(
        myvar.upper()) + '_regrid.nc'
    if ID >= 10: filename = dataPath + 'global-reanalysis-phys-001-009-ran-fr-glorys2-grid' + str(
        myvar.lower()) + '/REGRID/GLORYS2V1_ORCA025_' + str(year) + str(ID) + '15_R20110216_grid' + str(
        myvar.upper()) + '_regrid.nc'

    return filename


def getNORESMfilename(year, ID, myvar, dataPath):

    if (myvar=='grid'):
        filename = dataPath + 'GRID/NorESM.nc'
        #TODO: Fix this hardcoding of grid path
        #filename = "/work/users/trondk/REGSCEN/GRID/NorESM.nc"
    else:
        if myvar in ['iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']:
            if (ID < 10):
                filename = dataPath + 'ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h.'+str(year)+'-0'+str(ID)+'.nc'
            else:
                filename = dataPath + 'ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h.'+str(year)+'-'+str(ID)+'.nc'
        else:
            if (ID < 10):
                filename = dataPath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.'+str(year)+'-0'+str(ID)+'.nc'
            else:
                filename = dataPath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.'+str(year)+'-'+str(ID)+'.nc'

    print filename

    return filename

def getSODAMONTHLYfilename(year, ID, myvar, dataPath):

    if ID < 10: filename = dataPath + 'SODA_2.0.2_' + str(year) + '0' + str(ID) + '.cdf'
    if ID >= 10: filename = dataPath + 'SODA_2.0.2_' + str(year) + str(ID) + '.cdf'

    return filename

def getSODAfilename(year, ID, myvar, dataPath):

    file = "SODA_2.0.2_" + str(year) + "_" + str(ID) + ".cdf"

    return dataPath + file

def getWOAMONTHLYfilename(year, ID, myvar, dataPath):

    if myvar == "temperature":
        filename = dataPath + 'temperature_monthly_1deg.nc'
    elif myvar == "salinity":
        filename = dataPath + 'salinity_monthly_1deg.nc'
    else:
        print "Could not find any input files in folder: %s"%(datapath)

    return filename


def get3Ddata(grdROMS, grdMODEL, myvar, mytype, year, ID, varNames, dataPath):

    if myvar == 'temperature': varN = 0;
    if myvar == 'salinity': varN = 1;
    if myvar == 'uvel': varN = 3;
    if myvar == 'vvel': varN = 4;

    # The variable splitExtract is defined in IOsubset.py and depends on the orientation
    # and mytype of grid (-180-180 or 0-360). Assumes regular grid.
    if  grdMODEL.useESMF:
        if mytype == "SODA":
            filename = getSODAfilename(year, ID, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[varNames[varN]][0,:,:,:]

        if mytype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][:,:,:]

        if mytype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(year, ID, myvar, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][ID - 1, :,:,:]

        if mytype == "NORESM":
            cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))
            myunits = cdf.variables[str(varNames[varN])].units
            data = np.squeeze(cdf.variables[str(varNames[varN])][0,:,:,:])
            data = np.where(abs(data)>=32768 , grdROMS.fill_value, data)

        cdf.close()
    else:
        if grdMODEL.splitExtract is True:
            if mytype == "SODA":
                filename = getSODAfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[varNames[varN]][0, :,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[varNames[varN]][0, :,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if mytype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][:,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][:,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if mytype in ["WOAMONTHLY"]:

                filename = getWOAMONTHLYfilename(year, ID, myvar, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][ID - 1, :,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][ID - 1, :,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if mytype == "GLORYS2V1":
                if myvar == 'temperature':
                    cdf = Dataset(getGLORYS2V1filename(year, ID, 'T', dataPath))
                    myunits = cdf.variables[str(varNames[varN])].units

                if myvar == 'salinity': cdf = Dataset(getGLORYS2V1filename(year, ID, 'S', dataPath))
                if myvar == 'uvel': cdf = Dataset(getGLORYS2V1filename(year, ID, 'U', dataPath))
                if myvar == 'vvel': cdf = Dataset(getGLORYS2V1filename(year, ID, 'V', dataPath))

                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            if mytype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))
                myunits = cdf.variables[str(varNames[varN])].units
                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            cdf.close()

            data = np.concatenate((data1, data2), axis=2)

        else:
            if mytype == "SODA":
                filename = getSODAfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][0, :,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if mytype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][:,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if mytype == "WOAMONTHLY":
                filename = getWOAMONTHLYfilename(year, ID, myvar, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][ID, :,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if mytype == "GLORYS2V1":
                if myvar == 'temperature':
                    cdf = Dataset(getGLORYS2V1filename(year, ID, 'T', dataPath))

                    myunits = cdf.variables[str(varNames[varN])].units
                if myvar == 'salinity': cdf = Dataset(getGLORYS2V1filename(year, ID, 'S', dataPath))
                if myvar == 'uvel': cdf = Dataset(getGLORYS2V1filename(year, ID, 'U', dataPath))
                if myvar == 'vvel': cdf = Dataset(getGLORYS2V1filename(year, ID, 'V', dataPath))

                data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            if mytype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))
                myunits = cdf.variables[str(varNames[varN])].units

                data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
            cdf.close()

    if myvar == 'temperature' and mytype in ["GLORYS2V1", "NORESM"]:

        if myunits == "degree_Kelvin" or myunits == "K":
            if mytype in ["GLORYS2V1"]:
                data = np.where(data <= -32.767, grdROMS.fill_value, data)
            data = data - 273.15

            #  if time == 0 and myvar == myvars[0]:
            #      tmp = np.squeeze(data[0, :, :])
            #grdMODEL.mask = np.zeros(grdMODEL.lon.shape, dtype=np.float64)
            #grdMODEL.mask[:, :] = np.where(tmp == grdROMS.fill_value, 1, 0)

    if mytype == "GLORYS2V1":
        data = np.where(data <= -32.767, grdROMS.fill_value, data)
        data = np.ma.masked_where(data <= grdROMS.fill_value, data)


    if __debug__:
        print "Data range of %s just after extracting from netcdf file: %s - %s" % (str(varNames[varN]),
                                                                                    data.min(), data.max())

    return data

def get2Ddata(grdROMS, grdMODEL, myvar, mytype, year, ID, varNames, dataPath):

    indexROMS_SSH = (grdROMS.eta_rho, grdROMS.xi_rho)

    if mytype == "NORESM":
        if myvar == 'ageice': varN = 5;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'uice':   varN = 6;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'vice':   varN = 7;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'aice':   varN = 8;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':   varN = 9;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'snow_thick': varN = 10;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if myvar == 'ssh': varN = 2;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if  grdMODEL.useESMF:
        if mytype == "SODA":
            filename = getSODAfilename(year, ID, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[varNames[varN]][0,:,:]

        if mytype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][:,:]

        if mytype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(year, ID, myvar, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][ID - 1, :,:]

        if mytype == "NORESM":
            cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))
            #myunits = cdf.variables[str(varNames[varN])].units
            data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,:])

        cdf.close()
    else:
        if grdMODEL.splitExtract is True:
            if mytype == "SODA":
                filename = getSODAfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[varNames[varN]][0,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[varNames[varN]][0,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if mytype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if mytype == "GLORYS2V1":
                cdf = Dataset(getGLORYS2V1filename(year, ID, 'T', dataPath))

                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            if mytype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))

                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            cdf.close()
            data = np.concatenate((data1, data2), axis=1)

        else:
            if mytype == "SODA":
                filename = getSODAfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][0,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if mytype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if mytype == "GLORYS2V1":
                cdf = Dataset(getGLORYS2V1filename(year, ID, 'T', dataPath))

                data = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            if mytype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, ID, varNames[varN], dataPath))

                data = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            cdf.close()
            # data = np.where(abs(data) == grdMODEL.fill_value, grdROMS.fill_value, data)
            #data = np.ma.where(abs(data) > 1000, grdROMS.fill_value, data)

        if mytype == "GLORYS2V1":
            data = np.where(data <= -32.767, grdROMS.fill_value, data)
            data = np.ma.masked_where(data <= grdROMS.fill_value, data)

        if __debug__:
            print "Data range of %s just after extracting from netcdf file: %s - %s" % (str(varNames[varN]),
                                                                                        data.min(), data.max())

    return data


def convertMODEL2ROMS(years, IDS, climName, initName, dataPath, romsgridpath, myvars, show_progress, mytype, subset,
                      isClimatology, writeIce, useESMF):
    # First opening of input file is just for initialization of grid
    if mytype == 'SODA':
        fileNameIn = getSODAfilename(years[0], IDS[0], "salinity", dataPath)
    if mytype == 'SODAMONTHLY':
        fileNameIn = getSODAfilename(years[0], IDS[0], "salinity", dataPath)
    if mytype == 'NORESM':
        fileNameIn = getNORESMfilename(years[0], IDS[0], "grid", dataPath)
    if mytype == 'WOAMONTHLY':
        fileNameIn = getWOAMONTHLYfilename(years[0], IDS[0], "temperature", dataPath)
    if mytype == 'GLORYS2V1':
        fileNameIn = getGLORYS2V1filename(years[0], IDS[0], "S", dataPath)

    # First time in loop, get the essential old grid information
    # MODEL data already at Z-levels. No need to interpolate to fixed depths,
    # but we use the one we have
    grdMODEL = grd.grdClass(fileNameIn, mytype, useESMF)
    grdROMS = grd.grdClass(romsgridpath, "ROMS", useESMF)
    grdROMS.myvars = myvars
    if (useESMF):
        print "\nCreating the interpolation weights and indexes using ESMF:"
        print "  -> regridSrc2Dst at RHO points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_rho = ESMF.Field(grdROMS.esmfgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho, regrid_method=ESMF.RegridMethod.BILINEAR)
        print "  -> regridSrc2Dst at U points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_u = ESMF.Field(grdROMS.esmfgrid_u, "fieldDst_u", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_u = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_u, regrid_method=ESMF.RegridMethod.BILINEAR)
        print "  -> regridSrc2Dst at V points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_v = ESMF.Field(grdROMS.esmfgrid_v, "fieldDst_v", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_v = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_v, regrid_method=ESMF.RegridMethod.BILINEAR)

    # Now we want to subset the data to avoid storing more information than we need.
    # We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes
    IOsubset.findSubsetIndices(grdMODEL, min_lat=subset[0], max_lat=subset[1], min_lon=subset[2], max_lon=subset[3])

    print 'Initializing done'
    print '\n--------------------------'

    time = 0
    firstRun = True
    for year in years:

        for ID in IDS:

            if mytype == 'SODA':
                filename = getSODAfilename(year, ID, None, dataPath)
                varNames = ['TEMP', 'SALT', 'SSH', 'U', 'V']

            if mytype == 'SODAMONTHLY':
                filename = getSODAMONTHLYfilename(year, ID, None, dataPath)
                varNames = ['temp', 'salt', 'ssh', 'u', 'v']

            if mytype == 'GLORYS2V1':
                filename = getGLORYS2V1filename(year, ID, "S", dataPath)
                varNames = ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty']

            if mytype == 'WOAMONTHLY':
                filename = getWOAMONTHLYfilename(year, ID, "temperature", dataPath)
                varNames = ['t_an', 's_an']

            if mytype == 'NORESM':
                filename = getNORESMfilename(year, ID, "saln", dataPath)
                writeIce = True
                varNames = ['templvl','salnlvl','sealv', 'uvellvl', 'vvellvl','iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']

                myvars=['temperature','salinity', 'ssh', 'uvel', 'vvel','ageice','uice','vice','aice','hice','snow_thick']

                #varNames = ['iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']


            # Now open the input file and get the time
            cdf = Dataset(filename)

            getTime(cdf, grdROMS, grdMODEL, year, ID, time, mytype, firstRun)
            cdf.close()

            # Each MODEL file consist only of one time step. Get the subset data selected, and
            # store that time step in a new array:

            if firstRun is True:
                print "NOTICE!!: Make sure that these two arrays are in sequential order:"
                print "myvars:     %s" % (myvars)
                print "varnames: %s\n" % (varNames)
                firstRun = False

                # The first iteration we want to organize the subset indices we want to extract
                # from the input data to get the interpolation correct and to function fast
                IOsubset.organizeSplit(grdMODEL, grdROMS)

            for myvar in myvars:

                if myvar in ['temperature', 'salinity', 'uvel', 'vvel']:
                    data = get3Ddata(grdROMS, grdMODEL, myvar, mytype, year, ID, varNames, dataPath)

                if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice','snow_thick']:
                    data = get2Ddata(grdROMS, grdMODEL, myvar, mytype, year, ID, varNames, dataPath)

                # Take the input data and horizontally interpolate to your grid
                array1 = HorizontalInterpolation(myvar, grdROMS, grdMODEL, data, show_progress)

                if myvar in ['temperature', 'salinity']:
                    STdata = VerticalInterpolation(myvar, array1, array1, grdROMS, grdMODEL)
                    print "Data range of %s after interpolation: %3.3f to %3.3f" % (
                        myvar, STdata.min(), STdata.max())

                    for dd in xrange(len(STdata[:,0,0])):
                        STdata[dd,:,:] = np.where(grdROMS.mask_rho == 0, grdROMS.fill_value, STdata[dd,:,:])

                    STdata = np.where(abs(STdata) > 1000, grdROMS.fill_value, STdata)

                    IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, mytype, STdata)
                    if time == grdROMS.initTime and grdROMS.write_init is True:
                        IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, mytype, STdata)

                if myvar in ['ssh', 'ageice', 'aice', 'hice', 'snow_thick']:
                    SSHdata = array1[0, :, :]

                    SSHdata = np.where(grdROMS.mask_rho == 0, grdROMS.fill_value, SSHdata)
                    SSHdata = np.where(abs(SSHdata) > 100, grdROMS.fill_value, SSHdata)
                    SSHdata = np.where(abs(SSHdata) == 0, grdROMS.fill_value, SSHdata)

                    # Specific for ROMs. We set 0 where we should have fillvalue for ice otherwise ROMS blows up.
                    SSHdata = np.where(abs(SSHdata) == grdROMS.fill_value, 0, SSHdata)
                   # SSHdata = np.ma.masked_where(abs(SSHdata) > 100, SSHdata)

                    print "Data range of %s after interpolation: %3.3f to %3.3f" % (
                        myvar, SSHdata.min(), SSHdata.max())

                    IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, mytype, SSHdata)
                    if time == grdROMS.initTime:
                        IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, mytype, SSHdata)

                # The following are special routines used to calculate the u and v velocity
                # of ice based on the transport, which is divided by snow and ice thickenss
                # and then multiplied by grid size in dx or dy direction (opposite of transport).
                if myvar in ['uice', 'vice']:
                    SSHdata = array1[0, :, :]

                    if  myvar=="uice":mymask=grdROMS.mask_u
                    if  myvar=="vice":mymask=grdROMS.mask_v

                    SSHdata = np.where(mymask == 0, grdROMS.fill_value, SSHdata)
                    SSHdata = np.where(abs(SSHdata) > 100, grdROMS.fill_value, SSHdata)
                    SSHdata = np.where(abs(SSHdata) == 0, grdROMS.fill_value, SSHdata)
                    SSHdata = np.where(abs(SSHdata) == grdROMS.fill_value, 0, SSHdata)

                    #SSHdata = np.ma.masked_where(abs(SSHdata) > 1000, SSHdata)

                    print "Data range of %s after interpolation: %3.3f to %3.3f" % (myvar, SSHdata.min(), SSHdata.max())


                    IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, mytype, SSHdata)

                    if time == grdROMS.initTime:
                        if myvar == 'uice':
                            IOinitial.createInitFile(grdROMS, time, initName, 'uice', writeIce, mytype,  SSHdata)
                        if myvar == 'vice':
                            IOinitial.createInitFile(grdROMS, time, initName, 'vice', writeIce, mytype, SSHdata)

                if myvar == 'uvel':
                    array2 = array1

                if myvar == 'vvel':
                    indexROMS_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
                    indexROMS_VBAR = (grdROMS.eta_v, grdROMS.xi_v)
                    UBARdata = np.zeros((indexROMS_UBAR), dtype=np.float64)
                    VBARdata = np.zeros((indexROMS_VBAR), dtype=np.float64)

                    urot, vrot = rotate(grdROMS, grdMODEL, data, array2, array1)

                    u, v = interpolate2UV(grdROMS, grdMODEL, urot, vrot)

                    Udata, Vdata, UBARdata, VBARdata = VerticalInterpolation(myvar, u, v, grdROMS, grdMODEL)

                if myvar == 'vvel':
                    print "Data range of U after interpolation: %3.3f to %3.3f - V after scaling: %3.3f to %3.3f" % (
                        Udata.min(), Udata.max(), Vdata.min(), Vdata.max())

                    Udata = np.where(grdROMS.mask_u == 0, grdROMS.fill_value, Udata)
                    Udata = np.where(abs(Udata) > 1000, grdROMS.fill_value, Udata)
                    Vdata = np.where(grdROMS.mask_v == 0, grdROMS.fill_value, Vdata)
                    Vdata = np.where(abs(Vdata) > 1000, grdROMS.fill_value, Vdata)
                    UBARdata = np.where(grdROMS.mask_u == 0, grdROMS.fill_value, UBARdata)
                    UBARdata = np.where(abs(UBARdata) > 1000, grdROMS.fill_value, UBARdata)
                    VBARdata = np.where(grdROMS.mask_v == 0, grdROMS.fill_value, VBARdata)
                    VBARdata = np.where(abs(VBARdata) > 1000, grdROMS.fill_value, VBARdata)

                    IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, mytype, Udata, Vdata,
                                          UBARdata, VBARdata)
                    if time == grdROMS.initTime:
                        # We print time=initTime to init file so that we have values for ubar and vbar (not present at time=1)
                        IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, mytype, Udata, Vdata, UBARdata, VBARdata)

            time += 1
