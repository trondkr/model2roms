from __future__ import print_function

import datetime

import extrapolate as ex
import numpy as np

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2019, 3, 12)
__version__ = "1.5"
__status__ = "Development"


def laplacefilter(field, threshold, toxi, toeta):
    undef = 2.0e+35
    tx = 0.9 * undef
    critx = 0.01
    cor = 1.6
    mxs = 10

    field = np.where(abs(field) > threshold, undef, field)

    field = ex.extrapolate.fill(int(1), int(toxi),
                                int(1), int(toeta),
                                float(tx), float(critx), float(cor), float(mxs),
                                np.asarray(field, order='Fortran'),
                                int(toxi),
                                int(toeta))
    return field


def dohorinterpolationregulargrid(confM2R, mydata, myvar):
    if confM2R.show_progress is True:
        try:
            import progressbar
            widgets=['\rHorizontal interpolation:', progressbar.Percentage(), progressbar.Bar()]
            progress = progressbar.ProgressBar(confM2R.grdMODEL.nlevels, widgets=widgets).start()
        except ImportError:
            print("Could not find module progressbar")
            confM2R.show_progress=False
        pass

    indexROMS_Z_ST, toxi, toeta, mymask = setupIndexes(confM2R, myvar)
    array1 = np.zeros((indexROMS_Z_ST), dtype=np.float)

    # 2D or 3D interpolation
    depthlevels=confM2R.grdMODEL.nlevels
    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick','hs']:
        depthlevels=1
    
    for k in range(depthlevels):
        if confM2R.use_esmf:
            if depthlevels==1:
                indata=np.squeeze(mydata[:, :])
            else:
                indata=np.squeeze(mydata[k, :, :])

            # We interpolate to RHO fields for all variables and then we later interpolate RHO points to U and V points
            # But input data are read on U and V and RHO grids if they differ (as NorESM and GLORYS does).
            if myvar in ['uice','uvel']:
                confM2R.grdMODEL.fieldSrc_u.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_u(confM2R.grdMODEL.fieldSrc_u, confM2R.grdMODEL.fieldDst_rho)
            elif myvar in ['vice','vvel']:
                confM2R.grdMODEL.fieldSrc_v.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_v(confM2R.grdMODEL.fieldSrc_v, confM2R.grdMODEL.fieldDst_rho)
            else:
                confM2R.grdMODEL.fieldSrc_rho.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_rho(confM2R.grdMODEL.fieldSrc_rho, confM2R.grdMODEL.fieldDst_rho)
            
            # Since ESMF uses coordinates (x,y) we need to rotate and flip to get back to (y,x) order.
            field = np.fliplr(np.rot90(field.data, 3))
           
        if confM2R.use_filter:
            field = laplacefilter(field, 1000, toxi, toeta)
            field = field * mymask
            
        array1[k, :, :] = field

        if k in [39,34,17,2]:
            import plotData
            import matplotlib.pyplot as plt
            plotData.contourMap(confM2R.grdROMS, confM2R.grdROMS.lon_rho, confM2R.grdROMS.lat_rho, field, str(k)+'_withfilter', myvar)
            plotfilename="test_{}_wfilter.png".format(myvar)
            plt.savefig(plotfilename, dpi=150)
        # if __debug__:
        #      print "Data range after horisontal interpolation: ", field.min(), field.max()

        # if varname in ["hice","aice"]:
        #     import plotData
        #     plotData.contourMap(grdROMS,grdROMS.lon_rho,grdROMS.lat_rho, field, "surface", varname)

        if confM2R.show_progress is True:
            progress.update(k)
    if confM2R.show_progress is True:
        progress.finish()
    return array1

def setupIndexes(confM2R, myvar):
    if myvar in ["uice"]:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_u, confM2R.grdROMS.xi_u)
        toxi = confM2R.grdROMS.xi_u
        toeta = confM2R.grdROMS.eta_u
        mymask = confM2R.grdROMS.mask_u
    elif myvar in ["vice"]:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_v, confM2R.grdROMS.xi_v)
        toxi = confM2R.grdROMS.xi_v
        toeta = confM2R.grdROMS.eta_v
        mymask = confM2R.grdROMS.mask_v
    else:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_rho, confM2R.grdROMS.xi_rho)
        toxi = confM2R.grdROMS.xi_rho
        toeta = confM2R.grdROMS.eta_rho
        mymask = confM2R.grdROMS.mask_rho
    return indexROMS_Z_ST, toxi, toeta, mymask


def setupESMFInterpolationWeights(confM2R):
    if confM2R.use_esmf:
        print("=>Creating the interpolation weights and indexes using ESMF (this may take some time....):")

        print("  -> regridSrc2Dst at RHO points")
        confM2R.grdMODEL.fieldSrc_rho = ESMF.Field(confM2R.grdMODEL.esmfgrid, "fieldSrc", 
                                                    staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_rho = ESMF.Field(confM2R.grdROMS.esmfgrid, "fieldDst",
                                                   staggerloc=ESMF.StaggerLoc.CENTER)

        confM2R.grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_rho, 
                                                         confM2R.grdMODEL.fieldDst_rho,
                                                         regrid_method=ESMF.RegridMethod.BILINEAR,
                                                         unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("  -> regridSrc2Dst at U points to RHO")
        confM2R.grdMODEL.fieldSrc_u = ESMF.Field(confM2R.grdMODEL.esmfgrid_u, "fieldSrc", 
                                                       staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_u = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_u, 
                                                       confM2R.grdMODEL.fieldDst_rho,
                                                       regrid_method=ESMF.RegridMethod.BILINEAR,
                                                       unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("  -> regridSrc2Dst at V points to RHO")
        confM2R.grdMODEL.fieldSrc_v = ESMF.Field(confM2R.grdMODEL.esmfgrid_v, "fieldSrc", 
                                                       staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_v = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_v, 
                                                       confM2R.grdMODEL.fieldDst_rho,
                                                       regrid_method=ESMF.RegridMethod.BILINEAR,
                                                       unmapped_action=ESMF.UnmappedAction.IGNORE)