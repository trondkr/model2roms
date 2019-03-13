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
import forcingFilenames as fc

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass
__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2008, 8, 15)
__modified__ = datetime(2019, 3, 13)
__version__ = "1.8"
__status__ = "Development, modified on 15.08.2008,01.10.2009,07.01.2010, 15.07.2014, 01.12.2014, 07.08.2015, 08.02.2018, 04.03.2019, 13.03.2019"


def verticalinterpolation(myvar, array1, array2, grdROMS, grdMODEL):
    outINDEX_ST = (grdROMS.nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    outINDEX_U = (grdROMS.nlevels, grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_V = (grdROMS.nlevels, grdROMS.eta_v, grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v, grdROMS.xi_v)

    if myvar in ['salinity', 'temperature','O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
        print('\nStart vertical interpolation for %s (dimensions=%s x %s)' % (myvar, grdROMS.xi_rho, grdROMS.eta_rho))
        outdata = np.empty((outINDEX_ST), dtype=np.float, order='Fortran')
        
        outdata = interp.interpolation.dovertinter(np.asarray(outdata, order='F'),
                                                   np.asarray(array1, order='F'),
                                                   np.asarray(grdROMS.h, order='F'),
                                                   np.asarray(grdROMS.z_r, order='F'),
                                                   np.asarray(grdMODEL.z_r, order='F'),
                                                   grdROMS.nlevels,
                                                   grdMODEL.nlevels, 
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho))

        outdata = np.ma.masked_where(abs(outdata) > 1000, outdata)
        # The BCG has to be capped at 0
        if myvar in ['O3_c','O3_TA','N1_p','N3_p','N3_n','N5_s','O2_o']:
            outdata = np.ma.masked_where(abs(outdata) < 0, outdata)
       # import plotData
       # for k in range(len(grdROMS.nlevels)-1):
       #     plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, np.squeeze(outdata[k,:,:]),k, myvar)

        return outdata

    if myvar == 'vvel':
        print('\nStart vertical interpolation for uvel (dimensions=%s x %s)' % (grdROMS.xi_u, grdROMS.eta_u))
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

        print('\nStart vertical interpolation for vvel (dimensions=%s x %s)' % (grdROMS.xi_v, grdROMS.eta_v))
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
    Zu = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_u), int(grdROMS.xi_u)), np.float)
    Zv = np.zeros((int(grdMODEL.nlevels), int(grdROMS.eta_v), int(grdROMS.xi_v)), np.float)

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


def getTime(confM2R, year, month, day, ntime):
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    if confM2R.indatatype == 'SODA':
        filename = fc.getSODAfilename(confM2R, year, month, None)

    if confM2R.indatatype == 'SODA3':
        filename = fc.getSODA3filename(confM2R, year, month, None)

    if confM2R.indatatype == 'SODAMONTHLY':
        filename = fc.getSODAMONTHLYfilename(confM2R, year, month, None)

    if confM2R.indatatype == 'GLORYS':
        filename = fc.getGLORYSfilename(confM2R, year, month, "S")

    if confM2R.indatatype == 'WOAMONTHLY':
        filename = fc.getWOAMONTHLYfilename(confM2R, year, month, "temperature")

    if confM2R.indatatype == 'NORESM':
        filename = fc.getNORESMfilename(confM2R, year, month, "salnlvl")

    if confM2R.indatatype == 'NS8KM':
        filename = fc.getNS8KMfilename(confM2R, year, month, "salt")

    if confM2R.indatatype == 'NS8KMZ':
        filename, readFromOneFile = fc.getNS8KMZfilename(confM2R, year, month, "salt")

    # Now open the input file and get the time
    cdf = Dataset(filename)

    if confM2R.indatatype == 'NORESM':
        jdref = date2num(datetime(1948, 1, 1), units="days since 1948-01-01 00:00:00", calendar="standard")
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
        mydays = cdf.variables["time"][ntime]
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        # For NORESM we switch from in-units of 1800-01-01 to outunits of 1948-01-01
        currentdate = num2date(mydays, units=myunits, calendar=mycalendar)
        jd = date2num(currentdate,'days since 1948-01-01 00:00:00', calendar='noleap')

    confM2R.grdROMS.time = (jd - jdref)
    confM2R.grdROMS.reftime = jdref
    confM2R.grdROMS.timeunits = myunits
    cdf.close()
    print("-------------------------------")
    print('\nCurrent time of %s file : %s' % (confM2R.indatatype, currentdate))
    print("-------------------------------")

def get3ddata(confM2R, myvar, year, month, day, timecounter):
    varN=confM2R.globalvarnames.index(myvar)
   
    # The variable splitExtract is defined in IOsubset.py and depends on the orientation
    # and indatatype of grid (-180-180 or 0-360). Assumes regular grid.
    if confM2R.useesmf:
        filename = fc.getFilename(confM2R,year,month,confM2R.inputdatavarnames[varN])
        print(filename,confM2R.inputdatavarnames[varN])
        try:
            cdf = Dataset(filename)
        except:
            print("Unable to open input file {}".format(filename))
            return

        if confM2R.indatatype == "SODA":
            data = cdf.variables[confM2R.inputdatavarnames[varN]][0,:,:,:]

        if confM2R.indatatype == "SODA3":
            data = cdf.variables[confM2R.inputdatavarnames[varN]][month-1,:,:,:]

        if confM2R.indatatype == "SODAMONTHLY":
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][:, :, :]

        if confM2R.indatatype == "WOAMONTHLY":
            data = cdf.variables[str(confM2R.inputdatavarnames[varN])][month-1,:,:,:]

        if confM2R.indatatype == "NORESM":
            # For NorESM data - all data is in one big file so we need the timecounter to access correct data
            myunits = cdf.variables[str(confM2R.inputdatavarnames[varN])].units
            data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][timecounter,:,:,:])
            data = np.where(data.mask, confM2R.grdROMS.fillval, data)

        if confM2R.indatatype == "GLORYS":
            myunits = cdf.variables[str(confM2R.inputdatavarnames[varN])].units
            data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0,:,:,:])
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
        print("Data range of {} just after extracting from netcdf file: {:3.3f}-{:3.3f}".format(str(confM2R.inputdatavarnames[varN]),
                                                                                    float(data.min()), float(data.max())))
    return data


def get2ddata(confM2R, myvar, year, month, day, timecounter):
    
    varN=confM2R.globalvarnames.index(myvar)

    if confM2R.useesmf:
        
        if not confM2R.set2DvarsToZero:
            filename = fc.getFilename(confM2R,year,month,confM2R.inputdatavarnames[varN])
            try:
                cdf = Dataset(filename)
            except:
                print("Unable to open input file {}".format(filename))
                return

            if confM2R.indatatype == "SODA":
                data = cdf.variables[confM2R.inputdatavarnames[varN]][0,:,:]

            if confM2R.indatatype == "SODA3":
                if myvar == 'aice':
                    # We only extract the first thickness concentration. Need to fix this so all 5 classes can be extracted.
                    # http://www.atmos.umd.edu/~ocean/index_files/soda3_readme.htm
                    # hi: sea ice thickness [m ice]
                    # mi: sea ice mass [kg/m^2]
                    # hs: snow thickness [m snow]
                    # {cn1,cn2,cn3,cn4,cn5}: sea ice concentration [0:1] in five ice thickness classes
                    data = cdf.variables[confM2R.inputdatavarnames[varN]][int(month-1),0,:,:]
                else:
                    data = cdf.variables[confM2R.inputdatavarnames[varN]][int(month-1),:,:]

            if confM2R.indatatype == "SODAMONTHLY":
                data = cdf.variables[str(confM2R.inputdatavarnames[varN])][:, :]

            if confM2R.indatatype == "WOAMONTHLY":
                data = cdf.variables[str(confM2R.inputdatavarnames[varN])][month-1,:,:]

            if (confM2R.indatatype == "NORESM" and confM2R.set2DvarsToZero is False):
                # myunits = cdf.variables[str(grdROMS.varNames[varN])].units
                # For NORESM data are 12 months of data stored in ice files. Use ID as month indicator to get data.
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][timecounter,:,:])
                data = np.where(data.mask, confM2R.grdROMS.fillval, data)

            if confM2R.indatatype == "GLORYS":
                data = np.squeeze(cdf.variables[str(confM2R.inputdatavarnames[varN])][0,:,:])
                data = np.where(data.mask, confM2R.grdROMS.fillval, data)

            if not confM2R.set2DvarsToZero: cdf.close()

            if __debug__ and not confM2R.set2DvarsToZero:
                print("Data range of {} just after extracting from netcdf file: {:3.3f}-{:3.3f}".format(str(confM2R.inputdatavarnames[varN]),
                                                                                            float(data.min()), float(data.max())))
    if confM2R.set2DvarsToZero:
        return np.zeros((np.shape(confM2R.grdMODEL.lon)))
    return data


def convertMODEL2ROMS(confM2R):
    # First opening of input file is just for initialization of grid
    filenamein = fc.getFilename(confM2R,confM2R.start_year,confM2R.start_month,None)

    # Finalize creating the model grd object now that we know the filename for input data
    confM2R.grdMODEL.opennetcdf(filenamein)
    confM2R.grdMODEL.createobject(confM2R)
    confM2R.grdMODEL.getdims()
    # Create the ESMF weights used to do all of the horizontal interpolation
    interp2D.setupESMFInterpolationWeights(confM2R)
    
    # Now we want to subset the data to avoid storing more information than we need.
    # We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes
    if confM2R.subsetindata:
        IOsubset.findSubsetIndices(confM2R.grdMODEL, min_lat=confM2R.subset[0], max_lat=confM2R.subset[1],
                                   min_lon=confM2R.subset[2], max_lon=confM2R.subset[3])

    print('==> Initializing done')
    print('\n--------------------------')
    print('==> Starting loop over time')

    timecounter = 0
    firstrun = True

    for year in confM2R.years:
        months = datetimeFunctions.createlistofmonths(confM2R, year)

        for month in months:
            days = datetimeFunctions.createlistofdays(confM2R, year, month)

            for day in days:
                # Get the current date for given timestep 
                getTime(confM2R, year, month, day, timecounter)

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

                    if myvar in ['temperature','salinity','uvel','vvel','O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
                        data = get3ddata(confM2R, myvar, year, month, day, timecounter)

                    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
                        data = get2ddata(confM2R, myvar, year, month, day, timecounter)

                    # Take the input data and horizontally interpolate to your grid
                    array1 = interp2D.dohorinterpolationregulargrid(confM2R, data, myvar)

                    if myvar in ['temperature','salinity','O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
                        STdata = verticalinterpolation(myvar, array1, array1, confM2R.grdROMS, confM2R.grdMODEL)

                        for dd in range(len(STdata[:, 0, 0])):
                            STdata[dd, :, :] = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval,
                                                        STdata[dd, :, :])

                        STdata = np.where(abs(STdata) > 1000, confM2R.grdROMS.fillval, STdata)

                        IOwrite.writeclimfile(confM2R, timecounter, myvar, STdata)
                        if timecounter == confM2R.grdROMS.inittime and confM2R.grdROMS.write_init is True:
                            IOinitial.createinitfile(confM2R, timecounter, myvar, STdata)

                    if myvar in ['ssh','ageice','aice','hice','snow_thick']:
                        SSHdata = array1[0, :, :]

                        SSHdata = np.where(confM2R.grdROMS.mask_rho == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, confM2R.grdROMS.fillval, SSHdata)

                        # Specific for ROMs. We set 0 where we should have fillvalue for ice otherwise ROMS blows up.
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        IOwrite.writeclimfile(confM2R, timecounter,  myvar, SSHdata)

                        if timecounter == confM2R.grdROMS.inittime:
                            IOinitial.createinitfile(confM2R, timecounter,  myvar, SSHdata)

                    # The following are special routines used to calculate the u and v velocity
                    # of ice based on the transport, which is divided by snow and ice thickenss
                    # and then multiplied by grid size in dx or dy direction (opposite of transport).
                    if myvar in ['uice','vice']:
                        SSHdata = array1[0, :, :]

                        if myvar == "uice": mymask = confM2R.grdROMS.mask_u
                        if myvar == "vice": mymask = confM2R.grdROMS.mask_v

                        SSHdata = np.where(mymask == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, confM2R.grdROMS.fillval, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == confM2R.grdROMS.fillval, 0, SSHdata)

                        # SSHdata = np.ma.masked_where(abs(SSHdata) > 1000, SSHdata)

                     ##   print("Data range of {} after interpolation: {:3.3f}-{:3.3f}".format(
                      ##  myvar, float(SSHdata.min()), float(SSHdata.max())))

                        IOwrite.writeclimfile(confM2R, timecounter, myvar, SSHdata)

                        if timecounter == confM2R.grdROMS.inittime:
                            if myvar == 'uice':
                                IOinitial.createinitfile(confM2R, timecounter,  myvar, SSHdata)
                            if myvar == 'vice':
                                IOinitial.createinitfile(confM2R, timecounter,  myvar, SSHdata)

                    if myvar == 'uvel':
                        array2 = array1

                    if myvar == 'vvel':
                        urot, vrot = rotate(confM2R.grdROMS, confM2R.grdMODEL, data, array2, array1)
                        u, v = interpolate2uv(confM2R.grdROMS, confM2R.grdMODEL, urot, vrot)

                        Udata, Vdata, UBARdata, VBARdata = verticalinterpolation(myvar, u, v, confM2R.grdROMS,
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

                        IOwrite.writeclimfile(confM2R, timecounter, myvar,  Udata, Vdata, UBARdata, VBARdata)

                        if timecounter == confM2R.grdROMS.inittime:
                            IOinitial.createinitfile(confM2R, timecounter, myvar, Udata, Vdata, UBARdata, VBARdata)

                timecounter+=1
