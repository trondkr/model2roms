import os
import matplotlib

matplotlib.use('Agg')

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib


def contourMap(grdROMS, tlon, tlat, mydata, depthlevel, var):
    # TODO: Convert to cartopy
    fig = plt.figure(figsize=(12, 12), frameon=False)
    ax = plt.subplot(projection=ccrs.PlateCarree())
    steps = abs(np.max(mydata) - np.min(mydata)) / 20.
    levels = np.arange(np.min(mydata), np.max(mydata) + steps, steps)


    if var == 'temperature':
        levels = np.arange(-2.0, 20.0, 1)
    if var == 'salinity':
        levels = np.arange(0.0, 40.0, 0.3)
    if var == 'runoff':
        levels = np.arange(1e-4, 1e-8, 1e-10)
    if var == 'ssh':
        levels = np.arange(-1.5, 2.5, 0.25)
    if var in ['uvel', 'vvel', 'urot', 'vrot', 'ubar', 'vbar']:
        levels = np.arange(-2.0, 2.0, 0.1)
    if var == 'sodamask':
        levels = np.arange(0, 1, 1)
    if var == 'hice':
        levels = np.arange(0, 5, 0.1)
    if var == 'aice':
        levels = np.arange(0, 100, 5)
    if var == 'uice' or var == "vice":
        levels = np.arange(-1, 1, 0.1)
    if var == 'ageice':
        levels = np.arange(0, 12, 1)
    if var == 'snow_thick':
        levels = np.arange(0, 2, 0.1)

    ax.contourf(tlon, tlat, mydata, levels,
                transform=ccrs.PlateCarree())

    plt.title('Var:%s - depth:%s - time:%s' % (var, depthlevel, grdROMS.time))
    plotfile = 'figures/' + str(var) + '_depth_ESMF_' + str(depthlevel) + '_time_' + str(grdROMS.time) + '.png'
    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.savefig(plotfile)
    print("Saved figure: %s" % plotfile)
    plt.close()


def contourStationData(data, timedata, datedata, depthdata, stationName):
    x = timedata
    y = flipud(depthdata)
    X, Y = meshgrid(x, y)
    data = rot90(data)
    levels = np.arange(0, data.max(), 0.5)
    plt.figure(figsize=(10, 4))
    cs1 = contourf(X, Y, data, levels, cmap=cm.get_cmap('RdYlBu_r', len(levels) - 1), alpha=1.0, extend='both')
    cs2 = contour(X, Y, data, cs1.levels[::4], colors=('k',), hold='on', linewidths=(0.1,))
    colorbar(cs1, extend='both')

    sep = 73;
    s = 0;
    d = [];
    t = []

    for i in range(len(timedata)):
        if s == sep:
            d.append(datedata[i])
            t.append(timedata[i])
            s = 0

        s += 1
    locs, labels = plt.xticks(t, d)
    plt.setp(labels, 'rotation', 'vertical')

    plotfile = 'figures/station_' + str(stationName) + '.png'
    plt.savefig(plotfile)

    # plt.show()
