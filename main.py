import time, calendar
from netCDF4 import Dataset, datetime, date2num,num2date
import model2roms
import IOstation
import clim2bry
import decimateGrid
import grd
import numpy as np
import atmosForcing
import sys

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2017, 1, 3)
__version__ = "1.5"
__status__ = "Development"


def myhelp():
    """
    This program is run by typing: python main.py in the command window.
    """

def defineSubsetForIndata():
    # Subset the input data. The more you subset the less memory is needed for calculations
    # and the faster the process is performed. The subset is initially performed in IOsubset.py
    if gridtype == "NS8KM":
        minLat = 40; maxLat = 70; minLon = -30; maxLon = 40

    if gridtype == "A20":
        minLat = 30; maxLat = 90; minLon = -179; maxLon = 360

    if gridtype == "KINO":
        minLat = 30; maxLat = 70; minLon = -40; maxLon = 40

    if gridtype == "AA10KM":
        minLat = -70; maxLat = 89.5; minLon = -179; maxLon = 180

    subset = np.zeros(4); subset[0] = minLat; subset[1] = maxLat; subset[2] = minLon; subset[3] = maxLon
    return subset

def defineAbbreviation(gridtype):
    if gridtype == "NS8KM":
        return "nordsjoen_8km"

    if gridtype == "REGSCEN":
        return "regscen"
      
    if gridtype == "A20":
        return "a20"
      
    if gridtype == "KINO":
        return "kino"
      
    if gridtype == "AA10KM":
        return "aa10km"

    print "no abbreviation found"
    sys.exit()

def showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF, myformat):
    if isClimatology:
        print '\n=>Conversions run for climatological months'
    else:
        print '\n=>Conversions run from %s to year %s' % (start_year, end_year)
    print '==>The following variables will be converted:'
    for myvar in myvars:
        print '===> %s' % myvar
    if (useESMF):
        print "=>All horisontal interpolations will be done using ESMF-ESMPy (module ESMF)"
    print "=>Output files are written in format: %s"%(myformat)
    print '\n=>Output grid file is: %s' % romsgridpath
    

def formatDatesForOutputnames(start_year,end_year,start_month,end_month,start_day,end_day):
    # Format the date for use in output filenames
    startMonth = ("0%s"%(start_month) if start_month < 10 else "%s"%(start_month))
    startDay   = ("0%s"%(start_day) if start_day < 10 else "%s"%(start_day))
    endMonth   = ("0%s"%(end_month) if end_month < 10 else "%s"%(end_month))
    endDay     = ("0%s"%(end_day) if end_day < 10 else "%s"%(end_day))
    
    modelPeriod=str(start_year)+str(startMonth)+str(startDay) + '_to_' + str(end_year)+str(endMonth)+str(endDay)
   
    return modelPeriod

def defineOutputFilenames(abbreviation,start_year,end_year,start_month,end_month,start_day,end_day,indatatype):
    # Get string representation of start and end dates
    modelPeriod = formatDatesForOutputnames(start_year,end_year,start_month,end_month,start_day,end_day)
    
    # Name of output files for CLIM, BRY, and INIT files
    climName = abbreviation + '_clim_' + str(indatatype) + '_' + str(modelPeriod)+ '.nc'
    initName = abbreviation + '_init_' + str(indatatype) + '_' + str(modelPeriod)+ '.nc'
    bryName = abbreviation + '_bry_' + str(indatatype) + '_' + str(modelPeriod)+ '.nc'

    return climName, initName, bryName

def main():
    print '\n--------------------------\n'
    print 'Started ' + time.ctime(time.time())

    # EDIT ===================================================================
    # Set show_progress to "False" if you do not want to see the progress
    # indicator for horizontal interpolation.
    show_progress = True
    # Set compileAll to True if you want automatic re-compilation of all the
    # fortran files necessary to run model2roms. Options are "gfortran" or "ifort". Edit
    # compile.py to add other Fortran compilers.
    compileAll = False
    # Extract time-series of data for given longitude/latitude
    extractStations = False
     # Define a set of longitude/latitude positions with names to extract into
    # station files (using extractStations)
    if (extractStations):
      #  stationNames = ['NorthSea', 'Iceland', 'EastandWestGreenland', 'Lofoten', 'Georges Bank']
      #  lonlist = [2.4301, -22.6001, -47.0801, 13.3801, -67.2001]
      #  latlist = [54.5601, 63.7010, 60.4201, 67.5001, 41.6423]

        stationNames = ["Ytre Utsira","Indre Utsira","Lista"]
        latlist = [59.316667, 59.316667,58.016667]
        lonliost=[4.800000, 4.983333, 6.533333]

    # Create the bry, init, and clim files for a given grid and input data
    createOceanForcing = True
    # Create atmospheric forcing for the given grid
    createAtmosForcing = False # currently in beta stages and unavailable
    # Create a smaller resolution grid based on your original. Decimates every second for
    # each time run
    decimateGridfile = False
    # Write ice values to file (for Arctic regions)
    writeIce = True
    # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
    useESMF = True
    # Apply filter to smooth the 2D fields after interpolation (time consuming but enhances results)
    useFilter = True
    # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
    # Using NETCDF4 automatically turns on compression of files (ZLIB)
    myformat='NETCDF4'
    # Frequency of the input data: usually monthly 
    timeFrequencyOfInputData = "month" #, "month", "hour"

    # Subset input data. If you have global data you may want to seubset these to speed up reading. Make 
    # sure that your input data are cartesian (0-360 or -180:180, -90:90)
    subsetIndata = False
    if subsetIndata:
        subset = defineSubsetForIndata()

    # IN GRIDTYPES ------------------------------------------------------------------------------
    #  Define what grid type you wnat to interpolate from (input MODEL data)
    indatatype = 'SODA'
    indatatype = 'SODAMONTHLY'
    indatatype = 'WOAMONTHLY'
    indatatype = 'NORESM'
    indatatype = 'GLORYS'
    #indatatype = 'NS8KM'
    #indatatype = 'NS8KMZ'
    #indatatype = 'SODA3'

    # Define what grid type you wnat to interpolate from: Can be Z for SIGMA for ROMS
    # vertical coordinate system or ZLEVEL. also define the name of the dimensions in the input files.
    ingridtype = "SIGMA"
    ingridtype = "ZLEVEL"
    grdType  = 'regular'
    lonName = "nav_lon"
    latName = "nav_lat"
    depthName = "deptht"
    fill_value = -1.e+20

    # OUT GRIDTYPES ------------------------------------------------------------------------------
    # Define what grid type you wnat to interpolate to 
    outgrid  = "NS8KM"
    outgrid = "KINO"
    outgrid= "A20"
    outgridtype="ROMS"
    # Number of output depth levels
    Nlevels=40
    # define the grid stretching properties (leave default if uncertain what to pick)
    vstretching=4
    vtransform=2
    theta_s=7.0
    theta_b=0.1
    Tcline=250.0
    hc = 250
    
    # PATH TO FORCINGDATA --------------------------------------------------------------------
    # Define the path to the input data 
    if indatatype == 'SODA':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODA/"
    
    if indatatype == 'SODA3':
        modelpath = "/home/trondk/Projects/SODA3/SODA3/dsrs.atmos.umd.edu/DATA/soda3.3.1/REGRIDED/"
        modelpath = "/work/users/trondk/SODA3/dsrs.atmos.umd.edu/DATA/soda3.3.1/REGRIDED/"

    if indatatype == 'SODAMONTHLY':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODAMonthly/"
    
    if indatatype == 'GLORYS':
        modelpath = "/work/shared/imr/NS8KM/FORCING/GLORYS2V3/ftp.myocean.mercator-ocean.fr/Core/GLOBAL_REANALYSIS_PHYS_001_009/"
        
    if indatatype == 'NORESM':
        modelpath = "/Users/trondkr/Projects/RegScen/NRCP45AERCN_f19_g16_CLE_01/"
        if createAtmosForcing:
            atmospath = "/Users/trondkr/Projects/RegScen/model2roms/TESTFILES/"

    if indatatype == 'NS8KM':
        modelpath = "/Users/trondkr/Projects/is4dvar/grid2lonlat/RESULTS/"
    
    if indatatype == 'NS8KMZ':
        modelpath = "/Users/trondkr/Projects/NOWMAPS/delivery-COPERNICUS/"
        #modelpath = "/work/shared/imr/NS8KM/Z-LEVEL-ASSIMILATION2010-2013/"

    if indatatype == 'WOAMONTHLY':
        modelpath = "/Users/trondkr/Projects/is4dvar/createSSS/"


    # PATH TO GRID -----------------------------------------------------------------------------
    # Define the path to the grid file 
    if outgrid == "NS8KM":
        romsgridpath = "/work/users/trondk/NS8km/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v3.nc"

    if outgrid == "KINO":
        romsgridpath = "/work/shared/nn9297k/NORTHSEA/GRID/northsea_19022016_vf20.nc"
   
    if outgrid == "A20":
        romsgridpath="/home/trondk/Projects/A20/Grid/A20niva_grd_v1.nc"
        romsgridpath="/work/shared/nn9297k/A20/Grid/A20niva_grd_v1.nc"
        romsgridpath="/work/shared/nn9297k/FAABolous/WS4KM_grd.nc"
    
    if outgrid=="AA10KM":
        romsgridpath="/work-common/shared/imr/AA10KM/GRID/AA_10km_grid.nc"

    if indatatype == 'WOAMONTHLY': isClimatology = True
    else: isClimatology = False

    # DATE AND TIME DETAILS ---------------------------------------------------------
    # Define the period to create forcing for
    start_year  = 2013
    end_year    = 2015
    start_month = 3
    end_month   = 12
    start_day   = 15
    end_day     = 15

    if (int(calendar.monthrange(start_year, start_month)[1]) < start_day):
        start_day = int(calendar.monthrange(start_year, start_month)[1])

    if (int(calendar.monthrange(end_year, end_month)[1]) < end_day):
        end_day = int(calendar.monthrange(end_year, end_month)[1])

    startdate = datetime(start_year, start_month, start_day)
    enddate   = datetime(end_year, end_month, end_day)
    years = [start_year+year for year in xrange(end_year+1-start_year)]
   
    # DEFINE VARIABLE NAMES ---------------------------------------------------------
    # Define what and name of variables to include in the forcing files
    # -> myvars is the name model2roms uses to identify variables
    # -> varNames is the name of the variable found in the NetCDF input files
    if indatatype == 'SODA':
       myvars = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
       varNames = ['TEMP', 'SALT', 'SSH', 'U', 'V']

    if indatatype == 'SODA3':
       myvars = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel','aice','hice','snow_thick']
       varNames = ['temp', 'salt', 'ssh', 'u', 'v','cn', 'hi', 'hs']

    if indatatype == 'SODAMONTHLY':
        myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
        varNames = ['temp', 'salt', 'ssh', 'u', 'v']

    if indatatype ==  'NS8KM':
        myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
        varNames = ['votemper', 'vosaline', 'zeta', 'vozocrtx', 'vomecrty']
    
    if indatatype ==  'NS8KMZ':
        fileNameIn, readFromOneFile = model2roms.getNS8KMZfilename(startdate.year, startdate.month, startdate.day, "S", modelpath)
        
        myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']

        if (readFromOneFile):
             varNames = ['temp', 'salt', 'zeta', 'u_eastward', 'v_northward']
        else:
             varNames = ['votemper', 'vosaline', 'zeta', 'vozocrtx', 'vomecrty']

    if indatatype == 'GLORYS':
        if (writeIce):
            myvars   = ['temperature','salinity', 'ssh', 'uvel', 'vvel','uice','vice','aice','hice']
            varNames = ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty','iicevelu', 'iicevelv', 'ileadfra', 'iicethic']
        else:
            myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
            varNames = ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty']

    if indatatype == 'WOAMONTHLY':
        myvars   = ['temperature','salinity']
        varNames = ['t_an', 's_an']

    if indatatype == 'NORESM':
        myvars   = ['temperature','salinity', 'ssh', 'uvel', 'vvel','ageice','uice','vice','aice','hice','snow_thick']
        varNames = ['templvl','salnlvl','sealv', 'uvellvl', 'vvellvl','iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']


    # NO EDIT BELOW ====================================================================================================
    if compileAll is True:
        import compile; compile.compileAll("gfortran")

    if (createAtmosForcing or createOceanForcing):
        abbreviation = defineAbbreviation(outgrid)
        climName,initName,bryName = defineOutputFilenames(abbreviation,start_year,end_year,start_month,end_month,start_day,end_day,indatatype)

        if isClimatology is True:
            climName=abbreviation + '_' + str(indatatype) + '_climatology.nc'  

        showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF, myformat)

        # Create the grid object for the output grid
        grdROMS = grd.grdClass("ROMS", outgridtype, useESMF,'ocean', outgrid)
        grdROMS.Nlevels=Nlevels
        grdROMS.vstretching=vstretching
        grdROMS.vtransform=vtransform
        grdROMS.theta_s=theta_s
        grdROMS.theta_b=theta_b
        grdROMS.Tcline=Tcline
        grdROMS.hc = hc
        grdROMS.vars=myvars
        grdROMS.varNames=varNames
        grdROMS.lonName='lon_rho'
        grdROMS.latName='lat_rho'

        grdROMS.openNetCDF(romsgridpath)
        grdROMS.createObject()
        grdROMS.getDims()

        # Create the grid object for the input grid
        grdMODEL = grd.grdClass("FORCINGDATA", outgridtype, useESMF,'ocean', outgrid)
        grdMODEL.grdType  = grdType
        grdMODEL.lonName = lonName
        grdMODEL.latName = latName
        grdMODEL.depthName = depthName
        grdMODEL.fill_value = fill_value

        if (useESMF):
            import ESMF
            print "Sarting logfile for ESMF"
            manager = ESMF.Manager(logkind = ESMF.LogKind.MULTI, debug = True)

        if createOceanForcing:

            model2roms.convertMODEL2ROMS(grdMODEL,grdROMS,years, startdate, enddate, timeFrequencyOfInputData, climName, initName, modelpath, romsgridpath, show_progress,
                                             indatatype, outgridtype, isClimatology, writeIce, useESMF, useFilter, myformat, subsetIndata, outgrid, subset=None)

            clim2bry.writeBry(grdROMS, start_year, bryName, climName, writeIce, indatatype, myformat)

        if createAtmosForcing:
            atmosForcing.createAtmosFileUV(grdROMS,modelpath,atmospath,startdate,enddate,useESMF,
                myformat,abbreviation,indatatype,gridtype,show_progress)

    if decimateGridfile:
        decimateGrid.createGrid(grdROMS, "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015.nc", "/Users/trondkr/Projects/KINO/GRID/kino_1600m_18072015v2.nc", 2)

    if extractStations:
        print "Running in station mode and extracting pre-defined station locations"
        IOstation.getStationData(years, IDS, modelpath, latlist, lonlist, stationNames)

    print 'Finished ' + time.ctime(time.time())


if __name__ == "__main__":
    main()
