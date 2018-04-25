import numpy as np
import datetime
# import mpl_toolkits.basemap as mp
import extrapolate as ex

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2018, 4, 25)
__version__ = "1.5"
__status__ = "Development"


def laplaceFilter(field, threshold, toxi, toeta):
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


def doHorInterpolationRegularGrid(confM2R, mydata):
    if confM2R.showprogress is True:
        import progressbar
        # http://progressbar-2.readthedocs.org/en/latest/examples.html
        progress = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()],
                                           maxval=confM2R.grdMODEL.nlevels).start()
        # progress = progressbar.ProgressBar(widgets=[progressbar.BouncingBar(marker=progressbar.RotatingMarker(), fill_left=True)], maxval=grdMODEL.Nlevels).start()

    indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_rho, confM2R.grdROMS.xi_rho)
    array1 = np.zeros((indexROMS_Z_ST), dtype=np.float64)

    for k in range(confM2R.grdMODEL.nlevels):

        if confM2R.useesmf:
            confM2R.grdMODEL.fieldSrc.data[:, :] = np.flipud(np.rot90(np.squeeze(mydata[k, :, :])))
            # Get the actual regridded array
            field = confM2R.grdMODEL.regridSrc2Dst_rho(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_rho)

            # Since ESMF uses coordinates (x,y) we need to rotate and flip to get back to (y,x) order.
            field = np.fliplr(np.rot90(field.data, 3))

        if confM2R.usefilter:
            field = laplaceFilter(field, 1000, confM2R.grdROMS.xi_rho, confM2R.grdROMS.eta_rho)
        #  field=field*grdROMS.mask_rho

        array1[k, :, :] = field

        # if k in [34,17,2]:
        #     import plotData
        #     plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, field, str(k)+'_withfilter', myvar)
        # if __debug__:
        #      print "Data range after horisontal interpolation: ", field.min(), field.max()

        if confM2R.showprogress is True:
            progress.update(k)

    return array1


def doHorInterpolationSSHRegularGrid(confM2R, myvar, mydata):
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

    array1 = np.zeros((indexROMS_Z_ST), dtype=np.float64)

    if confM2R.useesmf:

        confM2R.grdMODEL.fieldSrc.data[:, :] = np.flipud(np.rot90(np.squeeze(mydata[:, :])))

        if myvar in ["uice"]:
            field = confM2R.grdMODEL.regridSrc2Dst_u(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_u)
        elif myvar in ["vice"]:
            field = confM2R.grdMODEL.regridSrc2Dst_v(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_v)
        else:
            field = confM2R.grdMODEL.regridSrc2Dst_rho(confM2R.grdMODEL.fieldSrc, confM2R.grdMODEL.fieldDst_rho)

        field = np.fliplr(np.rot90(field.data, 3))
        # if myvar in ["hice","aice"]:
        #     import plotData
        #     plotData.contourMap(grdROMS,grdROMS.lon_rho,grdROMS.lat_rho, field, "surface", myvar)

    # Smooth the output
    if confM2R.usefilter:
        field = laplaceFilter(field, 1000, toxi, toeta)
    field = field * mymask
    array1[0, :, :] = field

    #  import plotData
    #  plotData.contourMap(grdROMS, tolon, tolat, field, "34", myvar)

    return array1
