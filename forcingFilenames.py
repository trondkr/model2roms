import os
from netCDF4 import MFDataset, date2index, date2num, num2date
from datetime import datetime

# Main function called from model2roms
def getFilename(confM2R,year,month,day,defaultvar):
    if confM2R.ocean_indata_type == 'SODA':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getSODAfilename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == 'SODA3':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getSODA3filename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == 'SODA3_5DAY':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getSODA3_5DAYfilename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == 'SODAMONTHLY':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getSODAfilename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == 'NORESM':
        if defaultvar is None:defaultvar="grid"
        filenamein = getNORESMfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == 'WOAMONTHLY':
        if defaultvar is None:defaultvar="temperature"
        filenamein = getWOAMONTHLYfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == 'GLORYS':
        if defaultvar is None:defaultvar="S"
        filenamein = getGLORYSfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == 'GLORYS':
        if defaultvar is None:defaultvar="S"
        filenamein = getGLORYSfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == 'NS8KM':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getNS8KMfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == 'NS8KMZ':
        if defaultvar is None:defaultvar="salinity"
        filenamein, readFromOneFile = getNS8KMZfilename(confM2R, year, month, defaultvar)
    return filenamein
    
# private functions called from within module
def getGLORYSfilename(confM2R, year, month, myvar):
    # Month indicates month
    # myvar:S,T,U,V
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
        filename = confM2R.ocean_forcing_path + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' \
                   + str(myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(confM2R.start_year) + '0' + str(
            confM2R.start_month) \
                   + '15_' + str(production) + '_' + str(myvarPrefix) + '.nc'

    if confM2R.start_month >= 10:
        filename = confM2R.ocean_forcing_path + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' \
                   + str(myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(confM2R.start_year) + str(confM2R.start_month) \
                   + '15_' + str(production) + '_' + str(myvarPrefix) + '.nc'

    return filename


def getNORESMfilename(confM2R, year, month, myvar):
    if myvar == 'grid':
        filename = confM2R.ocean_forcing_path + 'grid_gx1v6.nc'
    else:
        if myvar in ['iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']:
            filename = confM2R.ocean_forcing_path + 'ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h.' + str(year) + '.nc'
        
        elif myvar in ['dissic','talk','po4','no3','si','o2']:
            filename = confM2R.ocean_forcing_path + "BCG_NRCP85BPRPEX_01.micom.hbgcmlvl.2006-2050.nc"
        
        elif myvar in ['templvl', 'salnlvl', 'sealv', 'uvellvl', 'vvellvl']:
            if myvar in ['salnlvl','templvl']:
                filename = confM2R.ocean_forcing_path + "TS_NRCP85BPRPEX_01.micom.2006-2100.nc"
            else:
                filename = confM2R.ocean_forcing_path + "VEL_NRCP85BPRPEX_01.micom.hmlvl.2006-2100.nc"
    return filename


def getNS8KMfilename(confM2R, year, month, myvar):
    if (confM2R.start_month < 10):
        mymonth = '0%s' % (confM2R.start_month)
    else:
        mymonth = '%s' % (confM2R.start_month)
    filename = confM2R.ocean_forcing_path + str(confM2R.start_year) + str(
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
        filename = confM2R.ocean_forcing_path + str(confM2R.start_year) + str(
            mymonth) + '15_mm-IMR-MODEL-ROMS-NWS-20160203-fv02.1.nc'

        return filename, readFromOneFile


def getSODAMONTHLYfilename(confM2R, year, month, myvar):
    if confM2R.start_month < 10:
        return confM2R.ocean_forcing_path + 'SODA_2.0.2_' + str(year) + '0' + str(month) + '.cdf'
    elif confM2R.start_month >= 10:
        return confM2R.ocean_forcing_path + 'SODA_2.0.2_' + str(year) + str(month) + '.cdf'


def getSODAfilename(confM2R, year, month, day, myvar):
    return confM2R.ocean_forcing_path + "SODA_2.0.2_" + str(year) + "_" + str(month) + ".cdf"

def getSODA3filename(confM2R, year, month, day, myvar):
    if (myvar in ['cn', 'hi', 'hs']):
        return confM2R.ocean_forcing_path + "soda3.3.1_mn_ice_reg_" + str(year) + ".nc"
    else:
        return confM2R.ocean_forcing_path + "soda3.3.1_mn_ocean_reg_" + str(year) + ".nc"

def getSODA3_5DAYfilename(confM2R, year, month, day, myvar):

    if len(confM2R.time_object)==0:
        mcdf = MFDataset(confM2R.ocean_forcing_path + "*.nc")
        confM2R.time_object = mcdf.variables["time"]
        
       # print("Loaded all timesteps: {}".format(confM2R.time_object[:]))
    index = date2index(datetime(year,month,day,0,0), confM2R.time_object, calendar=confM2R.time_object.calendar, select="nearest")
    seldate=num2date(confM2R.time_object[index], units=confM2R.time_object.units, calendar=confM2R.time_object.calendar)

    print("selected index {}".format(seldate))
    if (myvar in ['cn', 'hi', 'hs']):
        return '{}soda3.3.2_5dy_ocean_ice_{:04}_{:02}_{:02}.nc'.format(confM2R.ocean_forcing_path, seldate.year, seldate.month, seldate.day)
    else:
        return '{}soda3.3.2_5dy_ocean_reg_{:04}_{:02}_{:02}.nc'.format(confM2R.ocean_forcing_path, seldate.year, seldate.month, seldate.day)

def getERA5_1DAYfilename(confM2R, year, month, day, myvar):
    
    if len(confM2R.time_object)==0:
        mcdf = MFDataset(confM2R.atmospheric_forcing_path + "*.nc")
        confM2R.time_object = mcdf.variables["time"]
        
        print("Loaded all ERA5 timesteps: {}".format(confM2R.time_object[:]))
    index = date2index(datetime(year,month,day,0,0), confM2R.time_object, calendar=confM2R.time_object.calendar, select="nearest")
    seldate=num2date(confM2R.time_object[index], units=confM2R.time_object.units, calendar=confM2R.time_object.calendar)

    print("selected index {}".format(seldate))
 
    return '{}soda3.3.2_5dy_ocean_reg_{:04}_{:02}_{:02}.nc'.format(confM2R.atmospheric_forcing_path, seldate.year, seldate.month, seldate.day)



def getWOAMONTHLYfilename(confM2R, year, month, myvar):
    if myvar == "temperature":
        return confM2R.ocean_forcing_path + 'temperature_monthly_1deg.nc'
    elif myvar == "salinity":
        return confM2R.ocean_forcing_path + 'salinity_monthly_1deg.nc'
    else:
        print("Could not find any input files in folder: %s" % confM2R.ocean_forcing_path)