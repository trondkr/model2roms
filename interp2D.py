import numpy as np
import datetime
import mpl_toolkits.basemap as mp
import extrapolate as ex

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 4)
__modified__ = datetime.datetime(2014, 3, 11)
__version__  = "1.3"
__status__   = "Development"


def doHorInterpolationRegularGrid(myvar, grdROMS, grdMODEL, mydata, show_progress):

    if show_progress is True:
        import progressbar
        #http://progressbar-2.readthedocs.org/en/latest/examples.html
        progress = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], maxval=grdMODEL.Nlevels).start()
        #progress = progressbar.ProgressBar(widgets=[progressbar.BouncingBar(marker=progressbar.RotatingMarker(), fill_left=True)], maxval=grdMODEL.Nlevels).start()

    indexROMS_Z_ST = (grdMODEL.Nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    array1=np.zeros((indexROMS_Z_ST), dtype=np.float64)

    for k in xrange(grdMODEL.Nlevels):

        i0 = np.argmin(np.fabs(grdMODEL.lon[0,:]-180))

        mydataout = np.ma.zeros(np.squeeze(mydata[k,:,:]).shape,mydata.dtype)
        lonsout = np.ma.zeros((len(grdMODEL.lon[0,:])),grdMODEL.lon.dtype)

        lonsout[0:len(grdMODEL.lon[0,:])-i0] = grdMODEL.lon[0,i0:]-360

        lonsout[len(grdMODEL.lon[0,:])-i0:] = grdMODEL.lon[0,1:i0+1]

        mydataout[:,0:len(grdMODEL.lon[0,:])-i0]  = mydata[k,:,i0:]
        mydataout[:,len(grdMODEL.lon[0,:])-i0:] = mydata[k,:,1:i0+1]

        Zg = mp.interp(mydataout,lonsout,grdMODEL.lat[:,0],grdROMS.lon_rho,grdROMS.lat_rho,
                       checkbounds=False, masked=False, order=1)
      #  Zg=np.ma.masked_where(Zg==grdROMS.fill_value,Zg)

        undef=2.0e+35 #grdROMS.fill_value
        tx=0.9*undef
        critx=0.01
        cor=1.6
        mxs=100

        field = Zg
        field=np.where(abs(field) > 50 ,undef,field)

        field=ex.extrapolate.fill(int(1),int(grdROMS.xi_rho),
                                int(1),int(grdROMS.eta_rho),
                                float(tx), float(critx), float(cor), float(mxs),
                                np.asarray(field, order='Fortran'),
                                int(grdROMS.xi_rho),
                                int(grdROMS.eta_rho))


        field=field*grdROMS.mask_rho
        array1[k,:,:]=field

        #import plotmydata
        #if k==34:
        #    import plotData
        #    plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, field, k, myvar)
       # if __debug__:
       #     print "Data range after horisontal interpolation: ", field.min(), field.max()

        if show_progress is True:
           progress.update(k)

    return array1



def doHorInterpolationSSHRegularGrid(myvar,grdROMS,grdMODEL,mydata):

    if myvar in ["uice","sim"]:
        indexROMS_Z_ST = (grdMODEL.Nlevels,grdROMS.eta_u, grdROMS.xi_u)
        tolon = grdROMS.lon_u; tolat=grdROMS.lat_u
        toxi = grdROMS.xi_u; toeta=grdROMS.eta_u;
        mymask = grdROMS.mask_u
    elif myvar in ["vice","sim2"]:
        indexROMS_Z_ST = (grdMODEL.Nlevels,grdROMS.eta_v, grdROMS.xi_v)
        tolon = grdROMS.lon_v; tolat=grdROMS.lat_v
        toxi = grdROMS.xi_v; toeta=grdROMS.eta_v;
        mymask = grdROMS.mask_v
    else:
        indexROMS_Z_ST = (grdMODEL.Nlevels,grdROMS.eta_rho, grdROMS.xi_rho)
        tolon=grdROMS.lon_rho; tolat=grdROMS.lat_rho
        toxi=grdROMS.xi_rho; toeta=grdROMS.eta_rho;
        mymask = grdROMS.mask_rho

    array1=np.zeros((indexROMS_Z_ST),dtype=np.float64)

    i0 = np.argmin(np.fabs(grdMODEL.lon[0,:]-180))

    mydataout = np.zeros(np.squeeze(mydata[:,:]).shape,mydata.dtype)
    lonsout = np.zeros((len(grdMODEL.lon[0,:])),grdMODEL.lon.dtype)

    lonsout[0:len(grdMODEL.lon[0,:])-i0] = grdMODEL.lon[0,i0:]-360

    lonsout[len(grdMODEL.lon[0,:])-i0:] = grdMODEL.lon[0,1:i0+1]

    mydataout[:,0:len(grdMODEL.lon[0,:])-i0]  = mydata[:,i0:]
    mydataout[:,len(grdMODEL.lon[0,:])-i0:] = mydata[:,1:i0+1]

    Zg = mp.interp(mydataout,lonsout,grdMODEL.lat[:,0], tolon, tolat,
                   checkbounds=False, masked=False, order=1)

    undef=2.0e+35
    tx=0.9*undef
    critx=0.01
    cor=1.6
    mxs=100

    field = Zg
    field=np.where(abs(field)>100000,undef,field)

    field=ex.extrapolate.fill(int(1),int(toxi),
                           int(1),int(toeta),
                           float(tx), float(critx), float(cor), float(mxs),
                           np.asarray(field, order='Fortran'),
                           int(toxi),
                           int(toeta))

    field=field*mymask
    array1[0,:,:]=field

    #import plotData
    #plotData.contourMap(grdROMS, tolon, tolat, field, "34", myvar)

    return array1