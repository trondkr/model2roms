"""
This is CLM2BRY
"""
from datetime import datetime
from netCDF4 import Dataset, num2date
import numpy as np
import IOBry

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2009, 3, 2)
__modified__ = datetime(2019, 3, 12)
__version__ = "1.5"
__status__ = "Development"


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


def writebry(confM2R):
    # See myhelp function for definition of these variables:
    grdROMS = confM2R.grdROMS

    Lpo = grdROMS.Lp
    Mpo = grdROMS.Mp
    Lp = grdROMS.Lp - 1
    Mp = grdROMS.Mp - 1

    # Open the CLIM file
    clim = Dataset(confM2R.clim_name, 'r')
    # Generate the BRY netcdf4 file that we will use to fill in data
    IOBry.createBryFile(confM2R)
    # Now open the file we created
    f = Dataset(confM2R.bry_name, mode='a', format=confM2R.myformat, zlib=confM2R.myzlib)

    # Get the time from the clim file
    climtime = np.array(clim.variables["ocean_time"][:])
    ntimes = len(climtime)

    # For each time in CLIM file, save clips of boundary data to BRY file
    for itime in range(ntimes):

        temp = np.array(clim.variables["temp"][itime, :, :, :])
        salt = np.array(clim.variables["salt"][itime, :, :, :])
        ssh = np.array(clim.variables["zeta"][itime, :, :])
        u = np.array(clim.variables["u"][itime, :, :, :])
        v = np.array(clim.variables["v"][itime, :, :, :])
        ubar = np.array(clim.variables["ubar"][itime, :, :])
        vbar = np.array(clim.variables["vbar"][itime, :, :])

        if confM2R.writebcg:
            O3_c = np.array(clim.variables["O3_c"][itime, :, :, :])
            O3_TA = np.array(clim.variables["O3_TA"][itime, :, :, :])
            N1_p = np.array(clim.variables["N1_p"][itime, :, :, :])
            N3_n = np.array(clim.variables["N3_n"][itime, :, :, :])
            N5_s = np.array(clim.variables["N5_s"][itime, :, :, :])
            O2_o = np.array(clim.variables["O2_o"][itime, :, :, :])

        # Cut the boundary sections along Lpo, and Mpo to create the
        # North, South, West, and East boundary clips. Write each clip
        # to file for each time step

        # NOTE: For East and West only V current (v,vbar) have size Mp, others have size Mpo
        temp_west = np.squeeze(temp[:, :, 0])
        salt_west = np.squeeze(salt[:, :, 0])
        ssh_west = np.squeeze(ssh[:, 0])
        u_west = np.squeeze(u[:, :, 0])
        v_west = np.squeeze(v[:, :, 0])
        ubar_west = np.squeeze(ubar[:, 0])
        vbar_west = np.squeeze(vbar[:, 0])

        if confM2R.writebcg:
            O3_c_west = np.squeeze(O3_c[:, :, 0])
            O3_ta_west = np.squeeze(O3_TA[:, :, 0])
            N1_p_west = np.squeeze(N1_p[:, :, 0])
            N3_n_west = np.squeeze(N3_n[:, :, 0])
            N5_s_west = np.squeeze(N5_s[:, :, 0])
            O2_o_west = np.squeeze(O2_o[:, :, 0])

        temp_east = np.squeeze(temp[:, :, Lp])
        salt_east = np.squeeze(salt[:, :, Lp])
        ssh_east = np.squeeze(ssh[:, Lp])
        u_east = np.squeeze(u[:, :, Lp - 1])
        v_east = np.squeeze(v[:, :, Lp])
        ubar_east = np.squeeze(ubar[:, Lp - 1])
        vbar_east = np.squeeze(vbar[:, Lp])

        if confM2R.writebcg:
            O3_c_east = np.squeeze(O3_c[:, :, Lp])
            O3_ta_east = np.squeeze(O3_TA[:, :, Lp])
            N1_p_east = np.squeeze(N1_p[:, :, Lp])
            N3_n_east = np.squeeze(N3_n[:, :, Lp])
            N5_s_east = np.squeeze(N5_s[:, :, Lp])
            O2_o_east = np.squeeze(O2_o[:, :, Lp])

        # NOTE: For South and North only U current (u,ubar) have size Lp, others have size Lpo
        temp_south = np.squeeze(temp[:, 0, :])
        salt_south = np.squeeze(salt[:, 0, :])
        ssh_south = np.squeeze(ssh[0, :])
        u_south = np.squeeze(u[:, 0, :])
        v_south = np.squeeze(v[:, 0, :])
        ubar_south = np.squeeze(ubar[0, :])
        vbar_south = np.squeeze(vbar[0, :])

        if confM2R.writebcg:
            O3_c_south = np.squeeze(O3_c[:, 0, :])
            O3_ta_south = np.squeeze(O3_TA[:, 0, :])
            N1_p_south = np.squeeze(N1_p[:, 0, :])
            N3_n_south = np.squeeze(N3_n[:, 0, :])
            N5_s_south = np.squeeze(N5_s[:, 0, :])
            O2_o_south = np.squeeze(O2_o[:, 0, :])

        temp_north = np.squeeze(temp[:, Mp, :])
        salt_north = np.squeeze(salt[:, Mp, :])
        ssh_north = np.squeeze(ssh[Mp, :])
        u_north = np.squeeze(u[:, Mp, :])
        v_north = np.squeeze(v[:, Mp - 1, :])
        ubar_north = np.squeeze(ubar[Mp, :])
        vbar_north = np.squeeze(vbar[Mp - 1, :])

        if confM2R.writebcg:
            O3_c_north = np.squeeze(O3_c[:, Mp, :])
            O3_ta_north = np.squeeze(O3_TA[:, Mp, :])
            N1_p_north = np.squeeze(N1_p[:, Mp, :])
            N3_n_north = np.squeeze(N3_n[:, Mp, :])
            N5_s_north = np.squeeze(N5_s[:, Mp, :])
            O2_o_north = np.squeeze(O2_o[:, Mp, :])

        if confM2R.writeice:
            ageice = np.array(clim.variables["ageice"][itime, :, :])
            uice = np.array(clim.variables["uice"][itime, :, :])
            vice = np.array(clim.variables["vice"][itime, :, :])
            aice = np.array(clim.variables["aice"][itime, :, :])
            hice = np.array(clim.variables["hice"][itime, :, :])
            snow_thick = np.array(clim.variables["snow_thick"][itime, :, :])
            ti = np.array(clim.variables["ti"][itime, :, :])
            sfwat = np.array(clim.variables["sfwat"][itime, :, :])
            tisrf = np.array(clim.variables["tisrf"][itime, :, :])
            sig11 = np.array(clim.variables["sig11"][itime, :, :])
            sig12 = np.array(clim.variables["sig12"][itime, :, :])
            sig22 = np.array(clim.variables["sig22"][itime, :, :])

            ageice_west = np.squeeze(ageice[:, 0])
            ageice_east = np.squeeze(ageice[:, Lp])
            ageice_south = np.squeeze(ageice[0, :])
            ageice_north = np.squeeze(ageice[Mp, :])

            uice_west = np.squeeze(uice[:, 0])
            uice_east = np.squeeze(uice[:, Lp - 1])
            uice_south = np.squeeze(uice[0, :])
            uice_north = np.squeeze(uice[Mp, :])

            vice_west = np.squeeze(vice[:, 0])
            vice_east = np.squeeze(vice[:, Lp])
            vice_south = np.squeeze(vice[0, :])
            vice_north = np.squeeze(vice[Mp - 1, :])

            aice_west = np.squeeze(aice[:, 0])
            aice_east = np.squeeze(aice[:, Lp])
            aice_south = np.squeeze(aice[0, :])
            aice_north = np.squeeze(aice[Mp, :])

            hice_west = np.squeeze(hice[:, 0])
            hice_east = np.squeeze(hice[:, Lp])
            hice_south = np.squeeze(hice[0, :])
            hice_north = np.squeeze(hice[Mp, :])

            snow_thick_west = np.squeeze(snow_thick[:, 0])
            snow_thick_east = np.squeeze(snow_thick[:, Lp])
            snow_thick_south = np.squeeze(snow_thick[0, :])
            snow_thick_north = np.squeeze(snow_thick[Mp, :])

            ti_west = np.squeeze(ti[:, 0])
            ti_east = np.squeeze(ti[:, Lp])
            ti_south = np.squeeze(ti[0, :])
            ti_north = np.squeeze(ti[Mp, :])

            sfwat_west = np.squeeze(sfwat[:, 0])
            sfwat_east = np.squeeze(sfwat[:, Lp])
            sfwat_south = np.squeeze(sfwat[0, :])
            sfwat_north = np.squeeze(sfwat[Mp, :])

            tisrf_west = np.squeeze(tisrf[:, 0])
            tisrf_east = np.squeeze(tisrf[:, Lp])
            tisrf_south = np.squeeze(tisrf[0, :])
            tisrf_north = np.squeeze(tisrf[Mp, :])

            sig11_west = np.squeeze(sig11[:, 0])
            sig11_east = np.squeeze(sig11[:, Lp])
            sig11_south = np.squeeze(sig11[0, :])
            sig11_north = np.squeeze(sig11[Mp, :])

            sig12_west = np.squeeze(sig12[:, 0])
            sig12_east = np.squeeze(sig12[:, Lp])
            sig12_south = np.squeeze(sig12[0, :])
            sig12_north = np.squeeze(sig12[Mp, :])

            sig22_west = np.squeeze(sig22[:, 0])
            sig22_east = np.squeeze(sig22[:, Lp])
            sig22_south = np.squeeze(sig22[0, :])
            sig22_north = np.squeeze(sig22[Mp, :])

        # ------- Write time to file -------------------------------
        d = num2date(climtime[itime], units=clim.variables['ocean_time'].long_name,
                     calendar=clim.variables['ocean_time'].calendar)
        print('clim2bry.py => Appending data to file %s for time %s' % (confM2R.bry_name, d))

        f.variables['ocean_time'][itime] = climtime[itime]

        # ------- Write out western boundary variables ------------

        f.variables['temp_west'][itime, :, :] = temp_west
        f.variables['salt_west'][itime, :, :] = salt_west
        f.variables['zeta_west'][itime, :] = ssh_west

        f.variables['u_west'][itime, :, :] = u_west
        f.variables['v_west'][itime, :, :] = v_west
        f.variables['ubar_west'][itime, :] = ubar_west
        f.variables['vbar_west'][itime, :] = vbar_west

        if confM2R.writebcg:
            f.variables['O3_c_west'][itime, :, :] = O3_c_west
            f.variables['O3_TA_west'][itime, :, :] = O3_ta_west
            f.variables['N1_p_west'][itime, :, :] = N1_p_west
            f.variables['N3_n_west'][itime, :, :] = N3_n_west
            f.variables['N5_s_west'][itime, :, :] = N5_s_west
            f.variables['O2_o_west'][itime, :, :] = O2_o_west

        # ------- Write out eastern boundary variables ------------

        f.variables['temp_east'][itime, :, :] = temp_east
        f.variables['salt_east'][itime, :, :] = salt_east
        f.variables['zeta_east'][itime, :] = ssh_east
        f.variables['u_east'][itime, :, :] = u_east
        f.variables['v_east'][itime, :, :] = v_east
        f.variables['ubar_east'][itime, :] = ubar_east
        f.variables['vbar_east'][itime, :] = vbar_east

        if confM2R.writebcg:
            f.variables['O3_c_east'][itime, :, :] = O3_c_east
            f.variables['O3_TA_east'][itime, :, :] = O3_ta_east
            f.variables['N1_p_east'][itime, :, :] = N1_p_east
            f.variables['N3_n_east'][itime, :, :] = N3_n_east
            f.variables['N5_s_east'][itime, :, :] = N5_s_east
            f.variables['O2_o_east'][itime, :, :] = O2_o_east

        # ------- Write out southern boundary variables ------------

        f.variables['temp_south'][itime, :, :] = temp_south
        f.variables['salt_south'][itime, :, :] = salt_south
        f.variables['zeta_south'][itime, :] = ssh_south
        f.variables['u_south'][itime, :, :] = u_south
        f.variables['v_south'][itime, :, :] = v_south
        f.variables['ubar_south'][itime, :] = ubar_south
        f.variables['vbar_south'][itime, :] = vbar_south

        if confM2R.writebcg:
            f.variables['O3_c_south'][itime, :, :] = O3_c_south
            f.variables['O3_TA_south'][itime, :, :] = O3_ta_south
            f.variables['N1_p_south'][itime, :, :] = N1_p_south
            f.variables['N3_n_south'][itime, :, :] = N3_n_south
            f.variables['N5_s_south'][itime, :, :] = N5_s_south
            f.variables['O2_o_south'][itime, :, :] = O2_o_south

        # ------- Write out northern boundary variables ------------

        f.variables['temp_north'][itime, :, :] = temp_north
        f.variables['salt_north'][itime, :, :] = salt_north
        f.variables['zeta_north'][itime, :] = ssh_north
        f.variables['u_north'][itime, :, :] = u_north
        f.variables['v_north'][itime, :, :] = v_north
        f.variables['ubar_north'][itime, :] = ubar_north
        f.variables['vbar_north'][itime, :] = vbar_north

        if confM2R.writebcg:
            f.variables['O3_c_north'][itime, :, :] = O3_c_north
            f.variables['O3_TA_north'][itime, :, :] = O3_ta_north
            f.variables['N1_p_north'][itime, :, :] = N1_p_north
            f.variables['N3_n_north'][itime, :, :] = N3_n_north
            f.variables['N5_s_north'][itime, :, :] = N5_s_north
            f.variables['O2_o_north'][itime, :, :] = O2_o_north

        if confM2R.writeice:
            f.variables['ageice_west'][itime, :] = ageice_west
            f.variables['ageice_east'][itime, :] = ageice_east
            f.variables['ageice_south'][itime, :] = ageice_south
            f.variables['ageice_north'][itime, :] = ageice_north

            f.variables['uice_west'][itime, :] = uice_west
            f.variables['uice_east'][itime, :] = uice_east
            f.variables['uice_south'][itime, :] = uice_south
            f.variables['uice_north'][itime, :] = uice_north

            f.variables['vice_west'][itime, :] = vice_west
            f.variables['vice_east'][itime, :] = vice_east
            f.variables['vice_south'][itime, :] = vice_south
            f.variables['vice_north'][itime, :] = vice_north

            f.variables['aice_west'][itime, :] = aice_west
            f.variables['aice_east'][itime, :] = aice_east
            f.variables['aice_south'][itime, :] = aice_south
            f.variables['aice_north'][itime, :] = aice_north

            f.variables['hice_west'][itime, :] = hice_west
            f.variables['hice_east'][itime, :] = hice_east
            f.variables['hice_south'][itime, :] = hice_south
            f.variables['hice_north'][itime, :] = hice_north

            f.variables['snow_thick_west'][itime, :] = snow_thick_west
            f.variables['snow_thick_east'][itime, :] = snow_thick_east
            f.variables['snow_thick_south'][itime, :] = snow_thick_south
            f.variables['snow_thick_north'][itime, :] = snow_thick_north

            f.variables['ti_west'][itime, :] = ti_west
            f.variables['ti_east'][itime, :] = ti_east
            f.variables['ti_south'][itime, :] = ti_south
            f.variables['ti_north'][itime, :] = ti_north

            f.variables['sfwat_west'][itime, :] = sfwat_west
            f.variables['sfwat_east'][itime, :] = sfwat_east
            f.variables['sfwat_south'][itime, :] = sfwat_south
            f.variables['sfwat_north'][itime, :] = sfwat_north

            f.variables['tisrf_west'][itime, :] = tisrf_west
            f.variables['tisrf_east'][itime, :] = tisrf_east
            f.variables['tisrf_south'][itime, :] = tisrf_south
            f.variables['tisrf_north'][itime, :] = tisrf_north

            f.variables['sig11_west'][itime, :] = sig11_west
            f.variables['sig11_east'][itime, :] = sig11_east
            f.variables['sig11_south'][itime, :] = sig11_south
            f.variables['sig11_north'][itime, :] = sig11_north

            f.variables['sig12_west'][itime, :] = sig12_west
            f.variables['sig12_east'][itime, :] = sig12_east
            f.variables['sig12_south'][itime, :] = sig12_south
            f.variables['sig12_north'][itime, :] = sig12_north

            f.variables['sig22_west'][itime, :] = sig22_west
            f.variables['sig22_east'][itime, :] = sig22_east
            f.variables['sig22_south'][itime, :] = sig22_south
            f.variables['sig22_north'][itime, :] = sig22_north
    f.close()
