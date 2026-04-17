import sys
from datetime import datetime

import extrapolate as ex
import numpy as np
from netCDF4 import Dataset, num2date

import IOatmos
import grd
import forcingFilenames as fc

try:
    import ESMF
except:
    try:
        # The module name for ESMPy was changed in v8.4.0 from “ESMF” to “esmpy”
        import esmpy as ESMF
    except ImportError:
        raise ImportError("[atmosForcing]: Could not find module ESMF/esmpy. Required")
        sys.exit()

_author_   = 'Trond Kristiansen'
_email_    = 'me@trondkristiansen.com'
_created_  = datetime(2014, 12, 16)
_modified_ = datetime(2014, 12, 16)
_version_  = "0.2.0"
_status_   = "Development"


def help ():
    """
    This function creates atmospheric forcing files for ROMS

    def createAtmosFileUV(grdROMS, outfilename, output_format)

    To check the file for CF compliancy: http://titania.badc.rl.ac.uk/cgi-bin/cf-checker.pl?cfversion=1.0
    """


def laplaceFilter(field, threshold, toxi, toeta):
    undef=2.0e+35
    tx=0.9*undef
    critx=0.01
    cor=1.6
    mxs=10

    field=np.where(abs(field)>threshold,undef,field)

    field=ex.extrapolate.fill(int(1),int(toxi),
                                  int(1),int(toeta),
                                  float(tx), float(critx), float(cor), float(mxs),
                                  np.asarray(field, order='F'),
                                  int(toxi),
                                  int(toeta))
    return field


def getERA5Filename(confM2R):
    return confM2R.atmospheric_forcing_path + ''

def createAtmosFileUV(confM2R):
    
    if confM2R.show_progress is True:
        import progressbar
        progress = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], 
                                           maxval=len(confM2R.years)).start()

    # Prepare the objects for source and destination grids
    # ASSUMPTION: identical input grids for ocean and atmos variables
    if not confM2R.create_ocean_forcing:     # else already done in convert_MODEL2ROMS()
        # First opening of input file is just for initialization of grid
        filenamein = fc.get_filename(confM2R, confM2R.start_year, confM2R.start_month, confM2R.start_day, None)

        # Finalize creating the model grd object now that we know the filename for input data
        confM2R.grdMODEL.create_object(confM2R, filenamein)
        confM2R.grdMODEL.getdims()
    
    # Abbreviate grid objects
    grdROMS = confM2R.grdROMS
    grdMODEL = confM2R.grdMODEL
    
    # Create the outputfile
    outfilename= confM2R.abbreviation + '_windUV_' + str(confM2R.atmos_indata_type) + '_' \
                 + str(confM2R.startdate.year) + '_to_' + str(confM2R.enddate.year) + '.nc'
    IOatmos.createNetCDFFileUV(grdROMS, outfilename, confM2R.output_format, confM2R.atmos_indata_type)
    
    # Setup ESMF for interpolation (calculates weights)
    print("  -> regridSrc2Dst at RHO points")
    grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.fieldDst_rho = ESMF.Field(grdROMS.esmfgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho, 
                                             regrid_method=ESMF.RegridMethod.BILINEAR)

    if confM2R.atmos_indata_type == "NORESM":
        # compute wind speeds U,V from wind stress and interpolate
        # using variables TAUX, TAUY, U10
        year=2050; month=1; day=1           # only to get filename
        filename = fc.getNORESMfilename(year,month,day,"TAUX")
        cdf = Dataset(filename,"r")
        U10 = cdf.variables["U10"][:]       # 10 m wind speed
        TAUX = -(cdf.variables["TAUX"][:])  # zonal surface wind stress
        TAUY = -(cdf.variables["TAUY"][:])  # meridional surface wind stress

        magstr = np.sqrt(TAUX*TAUX + TAUY*TAUY)
        magstr = np.where(magstr < 1.e-8,1.e-8,magstr)

        windE = (TAUX/magstr)*U10
        windN = (TAUY/magstr)*U10
        
        time_in = cdf.variables["time"][:]
        time_calendar = cdf.variables['time'].calendar
        time_units = cdf.variables['time'].units

        scru = np.zeros((len(time_in),np.shape(grdROMS.lat_rho)[0],np.shape(grdROMS.lat_rho)[1]))
        scrv = np.zeros((len(time_in),np.shape(grdROMS.lat_rho)[0],np.shape(grdROMS.lat_rho)[1]))

        # Loop over each year and do the interpolations and write to file
        # Loop over each time-step in current file
        for t in range(len(time_in)):
            currentdate=num2date(time_in[t], units=time_units,calendar=time_calendar)
            print("Interpolating date: ",currentdate)

            # Eastward wind
            grdMODEL.fieldSrc[:,:]=np.flipud(np.rot90(np.squeeze(windE[t,:,:])))
            fieldE = grdMODEL.regridSrc2Dst_rho(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho)

            # Since ESMF uses coordinates (x,y) we need to rotate and flip to get back to (y,x) order.
            fieldE = np.fliplr(np.rot90(fieldE.data,3))
            fieldE = laplaceFilter(fieldE, 1000, grdROMS.xi_rho, grdROMS.eta_rho)
            fieldE = fieldE*grdROMS.mask_rho

            # Northward wind
            grdMODEL.fieldSrc[:,:]=np.flipud(np.rot90(np.squeeze(windN[t,:,:])))
            fieldN = grdMODEL.regridSrc2Dst_rho(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho)

            fieldN = np.fliplr(np.rot90(fieldN.data,3))
            fieldN = laplaceFilter(fieldN, 1000, grdROMS.xi_rho, grdROMS.eta_rho)
            fieldN = fieldN*grdROMS.mask_rho

            # Magnitude
            grdMODEL.fieldSrc[:,:]=np.flipud(np.rot90(np.squeeze(magstr[t,:,:])))
            magnitude = grdMODEL.regridSrc2Dst_rho(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho)

            magnitude = np.fliplr(np.rot90(magnitude.data,3))
            magnitude = laplaceFilter(magnitude, 1000, grdROMS.xi_rho, grdROMS.eta_rho)
            magnitude = magnitude*grdROMS.mask_rho

            import plotAtmos
            print("Interpolated range: ", np.min(magnitude), np.max(magnitude))
            print("Original range: ", np.min(magstr), np.max(magstr))

            grdROMS.time+=1

            print(np.shape(windE), np.shape(grdMODEL.lon), np.shape(grdMODEL.lat))
            plotAtmos.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, fieldE, fieldN, magnitude, 
                                 'wind','REGSCEN',currentdate)
            plotAtmos.contourMap(grdMODEL,
                                 grdMODEL.lon,
                                 grdMODEL.lat,
                                 np.squeeze(windE[t,:,:]),
                                 np.squeeze(windN[t,:,:]),
                                 np.squeeze(magstr[t,:,:]),
                                 'wind',confM2R.atmos_indata_type,currentdate)

            # Rotate to ROMS grid structure
            scru[t,:,:]=(fieldE*np.cos(grdROMS.angle)) + (fieldN*np.sin(grdROMS.angle))
            scrv[t,:,:]=(fieldN*np.cos(grdROMS.angle)) - (fieldE*np.sin(grdROMS.angle))
