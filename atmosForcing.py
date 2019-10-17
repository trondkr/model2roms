import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset, num2date
import numpy as np
import IOatmos
import grd
import extrapolate as ex

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF. Required")
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

    def createAtmosFileUV(grdROMS, outfilename, myformat)

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
                                  np.asarray(field, order='Fortran'),
                                  int(toxi),
                                  int(toeta))
    return field


def getERA5Filename(confM2R):
    return confM2R.atmosphericpath+''

def createAtmosFileUV(confM2R):
    
    if confM2R.showprogress is True:
        import progressbar
        progress = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], maxval=len(years)).start()
    
    # Create the objects for source and destination grids
   
    getERA5_1DAYfilename

    grdMODEL = grd.grdClass(nor, mytype, mytype, useESMF,'atmos')
    
    # Create the outputfile
    outfilename=  abbreviation + '_windUV_' + str(mytype) + '_' + str(startdate.year) + '_to_' + str(enddate.year) + '.nc'
    IOatmos.createNetCDFFileUV(grdROMS, outfilename, myformat, mytype)
    
    # Setup ESMF for interpolation (calculates weights)
    print("  -> regridSrc2Dst at RHO points")
    grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.fieldDst_rho = ESMF.Field(grdROMS.esmfgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho, regrid_method=ESMF.RegridMethod.BILINEAR)

    # Loop over each year and do the interpolations and write to file
    year=2050; month=1; day=1
    if mytype == "NORESM":
        
        filename = getNORESMfilename(year,month,day,"TAUX",atmospath)
        cdf = Dataset(filename,"r")
        U10 = cdf.variables["U10"][:]
        TAUX = -(cdf.variables["TAUX"][:])
        TAUY = -(cdf.variables["TAUY"][:])

        magstr = np.sqrt(TAUX*TAUX + TAUY*TAUY)
        magstr = np.where(magstr < 1.e-8,1.e-8,magstr)

        windE = (TAUX/magstr)*U10
        windN = (TAUY/magstr)*U10
        
        time_in = cdf.variables["time"][:]
        time_calendar = cdf.variables['time'].calendar
        time_units = cdf.variables['time'].units
       
        scru = np.zeros((len(time_in),np.shape(grdROMS.lat_rho)[0],np.shape(grdROMS.lat_rho)[1]))
        scrv = np.zeros((len(time_in),np.shape(grdROMS.lat_rho)[0],np.shape(grdROMS.lat_rho)[1]))
  
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
                'wind','NORESM',currentdate)

            # Rotate to ROMS grid structure
            scru[t,:,:]=(fieldE*np.cos(grdROMS.angle)) + (fieldN*np.sin(grdROMS.angle))
            scrv[t,:,:]=(fieldN*np.cos(grdROMS.angle)) - (fieldE*np.sin(grdROMS.angle))
         

