
"""
This is CLM2BRY
"""
from datetime import datetime
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import IOBry

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 3, 2)
__modified__ = datetime(2014, 4, 7)
__version__  = "1.2"
__status__   = "Development"

def myhelp():
    """
    This script generates boundary (BRY) files from the climatology (CLIM) files. The
    climatology files are created using createForcing option in main.py of soda2roms package.

    Since the variables have different lengths in eta and xi directions, the
    clips of data along the East, WEst, and North and South transects will differ in size. The correct
    sizes are defined below, where No is the number of vertical levels (length of s_rho):

    Define the sizes of the clips along xi :
    North and South = salt(No,Lpo),
                      temp(No,Lpo),
                      v(No,Lpo),
                      vbar(1,Lpo),
                      zeta(1,Lpo),
                      u(No,Lp)
                      ubar(1,Lp)

    Define the sizes of the clips along eta :
    East and West   = salt(No,Mpo),
                      temp(No,Mpo),
                      v(No,Mp),
                      vbar(1,Mp),
                      zeta(1,Mpo),
                      u(No,Mpo)
                      ubar(1,Mpo)"""

def writeBry(grdROMS, year, bryName, climName, writeIce, mytype):

    myzlib = True
    myformat='NETCDF4'

    # See myhelp function for definition of these variables:

    Lpo = grdROMS.Lp
    Mpo = grdROMS.Mp
    Lp=grdROMS.Lp-1
    Mp=grdROMS.Mp-1

    # Open the CLIM file
    clim    = Dataset(climName,'r')
    # Generate the BRY netcdf4 file that we will use to fill in data
    IOBry.createBryFile(grdROMS, bryName, writeIce, mytype)
    # Now open the file we created
    f = Dataset(bryName, mode='a', format=myformat, zlib=myzlib)

    # Get the time from the clim file
    climtime = np.array(clim.variables["ocean_time"][:])
    ntimes = len(climtime)

    # For each time in CLIM file, save clips of boundary data to BRY file
    for itime in range(ntimes):

        temp        = np.array(clim.variables["temp"][itime,:,:,:])
        salt        = np.array(clim.variables["salt"][itime, :,:,:])
        ssh         = np.array(clim.variables["zeta"][itime,:,:])
        u           = np.array(clim.variables["u"][itime,:,:,:])
        v           = np.array(clim.variables["v"][itime,:,:,:])
        ubar        = np.array(clim.variables["ubar"][itime,:,:])
        vbar        = np.array(clim.variables["vbar"][itime,:,:])

        # Cut the boundary sections along Lpo, and Mpo to create the
        # North, South, West, and East boundary clips. Write each clip
        # to file for each time step

        # NOTE: For East and West only V current (v,vbar) have size Mp, others have size Mpo
        temp_west = np.squeeze(temp[:,:,0])
        salt_west = np.squeeze(salt[:,:,0])
        ssh_west  = np.squeeze(ssh[:,0])
        u_west    = np.squeeze(u[:,:,0])
        v_west    = np.squeeze(v[:,:,0])
        ubar_west = np.squeeze(ubar[:,0])
        vbar_west = np.squeeze(vbar[:,0])

        temp_east = np.squeeze(temp[:,:,Lp])
        salt_east = np.squeeze(salt[:,:,Lp])
        ssh_east  = np.squeeze(ssh[:,Lp])
        u_east    = np.squeeze(u[:,:,Lp-1])
        v_east    = np.squeeze(v[:,:,Lp])
        ubar_east = np.squeeze(ubar[:,Lp-1])
        vbar_east = np.squeeze(vbar[:,Lp])

        # NOTE: For South and North only U current (u,ubar) have size Lp, others have size Lpo
        temp_south = np.squeeze(temp[:,0,:])
        salt_south = np.squeeze(salt[:,0,:])
        ssh_south  = np.squeeze(ssh[0,:])
        u_south    = np.squeeze(u[:,0,:])
        v_south    = np.squeeze(v[:,0,:])
        ubar_south = np.squeeze(ubar[0,:])
        vbar_south = np.squeeze(vbar[0,:])

        temp_north = np.squeeze(temp[:,Mp,:])
        salt_north = np.squeeze(salt[:,Mp,:])
        ssh_north  = np.squeeze(ssh[Mp,:])
        u_north    = np.squeeze(u[:,Mp,:])
        v_north    = np.squeeze(v[:,Mp-1,:])
        ubar_north = np.squeeze(ubar[Mp,:])
        vbar_north = np.squeeze(vbar[Mp-1,:])

        if writeIce:
            ageice     = np.array(clim.variables["ageice"][itime,:,:])
            uice       = np.array(clim.variables["uice"][itime,:,:])
            vice       = np.array(clim.variables["vice"][itime,:,:])
            aice       = np.array(clim.variables["aice"][itime,:,:])
            hice       = np.array(clim.variables["hice"][itime,:,:])
            snow_thick = np.array(clim.variables["snow_thick"][itime,:,:])
            ti         = np.array(clim.variables["ti"][itime,:,:])
            sfwat      = np.array(clim.variables["sfwat"][itime,:,:])
            tisrf      = np.array(clim.variables["tisrf"][itime,:,:])
            sig11      = np.array(clim.variables["sig11"][itime,:,:])
            sig12      = np.array(clim.variables["sig12"][itime,:,:])
            sig22      = np.array(clim.variables["sig22"][itime,:,:])

            ageice_west  = np.squeeze(ageice[:,0])
            ageice_east  = np.squeeze(ageice[:,Lp])
            ageice_south = np.squeeze(ageice[0,:])
            ageice_north = np.squeeze(ageice[Mp,:])

            uice_west  = np.squeeze(uice[:,0])
            uice_east  = np.squeeze(uice[:,Lp-1])
            uice_south = np.squeeze(uice[0,:])
            uice_north = np.squeeze(uice[Mp,:])

            vice_west  = np.squeeze(vice[:,0])
            vice_east  = np.squeeze(vice[:,Lp])
            vice_south = np.squeeze(vice[0,:])
            vice_north = np.squeeze(vice[Mp-1,:])

            aice_west  = np.squeeze(aice[:,0])
            aice_east  = np.squeeze(aice[:,Lp])
            aice_south = np.squeeze(aice[0,:])
            aice_north = np.squeeze(aice[Mp,:])

            hice_west  = np.squeeze(hice[:,0])
            hice_east  = np.squeeze(hice[:,Lp])
            hice_south = np.squeeze(hice[0,:])
            hice_north = np.squeeze(hice[Mp,:])

            snow_thick_west  = np.squeeze(snow_thick[:,0])
            snow_thick_east  = np.squeeze(snow_thick[:,Lp])
            snow_thick_south = np.squeeze(snow_thick[0,:])
            snow_thick_north = np.squeeze(snow_thick[Mp,:])

            ti_west  = np.squeeze(ti[:,0])
            ti_east  = np.squeeze(ti[:,Lp])
            ti_south = np.squeeze(ti[0,:])
            ti_north = np.squeeze(ti[Mp,:])

            sfwat_west  = np.squeeze(sfwat[:,0])
            sfwat_east  = np.squeeze(sfwat[:,Lp])
            sfwat_south = np.squeeze(sfwat[0,:])
            sfwat_north = np.squeeze(sfwat[Mp,:])

            tisrf_west  = np.squeeze(tisrf[:,0])
            tisrf_east  = np.squeeze(tisrf[:,Lp])
            tisrf_south = np.squeeze(tisrf[0,:])
            tisrf_north = np.squeeze(tisrf[Mp,:])

            sig11_west  = np.squeeze(sig11[:,0])
            sig11_east  = np.squeeze(sig11[:,Lp])
            sig11_south = np.squeeze(sig11[0,:])
            sig11_north = np.squeeze(sig11[Mp,:])

            sig12_west  = np.squeeze(sig12[:,0])
            sig12_east  = np.squeeze(sig12[:,Lp])
            sig12_south = np.squeeze(sig12[0,:])
            sig12_north = np.squeeze(sig12[Mp,:])

            sig22_west  = np.squeeze(sig22[:,0])
            sig22_east  = np.squeeze(sig22[:,Lp])
            sig22_south = np.squeeze(sig22[0,:])
            sig22_north = np.squeeze(sig22[Mp,:])

        # ------- Write time to file -------------------------------
        d= num2date(climtime[itime],units=clim.variables['ocean_time'].long_name,calendar=clim.variables['ocean_time'].calendar)
        print 'clim2bry.py => Appending data to file %s for time %s'%(bryName,d)

        f.variables['ocean_time'][itime]      = climtime[itime]

        # ------- Write out western boundary variables ------------

        f.variables['temp_west'][itime,:,:] = temp_west
        f.variables['salt_west'][itime,:,:] = salt_west
        f.variables['zeta_west'][itime,:]   = ssh_west

        f.variables['u_west'][itime,:,:]    = u_west
        f.variables['v_west'][itime,:,:]    = v_west
        f.variables['ubar_west'][itime,:]   = ubar_west
        f.variables['vbar_west'][itime,:]   = vbar_west

        # ------- Write out eastern boundary variables ------------

        f.variables['temp_east'][itime,:,:] = temp_east
        f.variables['salt_east'][itime,:,:] = salt_east
        f.variables['zeta_east'][itime,:]   = ssh_east
        f.variables['u_east'][itime,:,:]    = u_east
        f.variables['v_east'][itime,:,:]    = v_east
        f.variables['ubar_east'][itime,:]   = ubar_east
        f.variables['vbar_east'][itime,:]   = vbar_east

        # ------- Write out southern boundary variables ------------

        f.variables['temp_south'][itime,:,:] = temp_south
        f.variables['salt_south'][itime,:,:] = salt_south
        f.variables['zeta_south'][itime,:]   = ssh_south
        f.variables['u_south'][itime,:,:]    = u_south
        f.variables['v_south'][itime,:,:]    = v_south
        f.variables['ubar_south'][itime,:]   = ubar_south
        f.variables['vbar_south'][itime,:]   = vbar_south

        # ------- Write out northern boundary variables ------------

        f.variables['temp_north'][itime,:,:] = temp_north
        f.variables['salt_north'][itime,:,:] = salt_north
        f.variables['zeta_north'][itime,:]   = ssh_north
        f.variables['u_north'][itime,:,:]    = u_north
        f.variables['v_north'][itime,:,:]    = v_north
        f.variables['ubar_north'][itime,:]   = ubar_north
        f.variables['vbar_north'][itime,:]   = vbar_north

        if writeIce:
            f.variables['ageice_west'][itime,:] = ageice_west
            f.variables['ageice_east'][itime,:] = ageice_east
            f.variables['ageice_south'][itime,:] = ageice_south
            f.variables['ageice_north'][itime,:] = ageice_north

            f.variables['uice_west'][itime,:] = uice_west
            f.variables['uice_east'][itime,:] = uice_east
            f.variables['uice_south'][itime,:] = uice_south
            f.variables['uice_north'][itime,:] = uice_north

            f.variables['vice_west'][itime,:] = vice_west
            f.variables['vice_east'][itime,:] = vice_east
            f.variables['vice_south'][itime,:] = vice_south
            f.variables['vice_north'][itime,:] = vice_north

            f.variables['aice_west'][itime,:] = aice_west
            f.variables['aice_east'][itime,:] = aice_east
            f.variables['aice_south'][itime,:] = aice_south
            f.variables['aice_north'][itime,:] = aice_north

            f.variables['hice_west'][itime,:] = hice_west
            f.variables['hice_east'][itime,:] = hice_east
            f.variables['hice_south'][itime,:] = hice_south
            f.variables['hice_north'][itime,:] = hice_north

            f.variables['snow_thick_west'][itime,:] = snow_thick_west
            f.variables['snow_thick_east'][itime,:] = snow_thick_east
            f.variables['snow_thick_south'][itime,:] = snow_thick_south
            f.variables['snow_thick_north'][itime,:] = snow_thick_north

            f.variables['ti_west'][itime,:] = ti_west
            f.variables['ti_east'][itime,:] = ti_east
            f.variables['ti_south'][itime,:] = ti_south
            f.variables['ti_north'][itime,:] = ti_north

            f.variables['sfwat_west'][itime,:] = sfwat_west
            f.variables['sfwat_east'][itime,:] = sfwat_east
            f.variables['sfwat_south'][itime,:] = sfwat_south
            f.variables['sfwat_north'][itime,:] = sfwat_north

            f.variables['tisrf_west'][itime,:] = tisrf_west
            f.variables['tisrf_east'][itime,:] = tisrf_east
            f.variables['tisrf_south'][itime,:] = tisrf_south
            f.variables['tisrf_north'][itime,:] = tisrf_north

            f.variables['sig11_west'][itime,:] = sig11_west
            f.variables['sig11_east'][itime,:] = sig11_east
            f.variables['sig11_south'][itime,:] = sig11_south
            f.variables['sig11_north'][itime,:] = sig11_north

            f.variables['sig12_west'][itime,:] = sig12_west
            f.variables['sig12_east'][itime,:] = sig12_east
            f.variables['sig12_south'][itime,:] = sig12_south
            f.variables['sig12_north'][itime,:] = sig12_north

            f.variables['sig22_west'][itime,:] = sig22_west
            f.variables['sig22_east'][itime,:] = sig22_east
            f.variables['sig22_south'][itime,:] = sig22_south
            f.variables['sig22_north'][itime,:] = sig22_north
    f.close()
