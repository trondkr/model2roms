from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile
from pylab import *

def contourMap(grdROMS,grdSODA,data,depthlevel,var):
    
    
    tlat      = np.array(grdROMS.lat_rho)
    tlon      = np.array(grdROMS.lon_rho)
    
    if var=='vrot':
        tlat      = np.array(grdROMS.lat_v)
        tlon      = np.array(grdROMS.lon_v)
    if var=='urot':
        tlat      = np.array(grdROMS.lat_u)
        tlon      = np.array(grdROMS.lon_u)
    
    temp      = data
    
    
    plt.figure(figsize=(12,12))

    if grdROMS.grdName=='GOM':
    # Plot the Georges Bank region
        map = Basemap(llcrnrlon=-80.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
                      resolution='i',projection='tmerc',lon_0=-70,lat_0=0,area_thresh=10.)
        levels = np.arange(0.0, 25.0, 0.5)
        
    if grdROMS.grdName=='Nordic':
    # Plot the Nordic region (Greenland, Nordic Seas, and the Barents Sea)
        map = Basemap(lon_0=25,boundinglat=50,
                      resolution='l',area_thresh=500.,projection='npstere')
        levels = np.arange(-2.0, 13.0, 0.1)
        
    if grdROMS.grdName=='NA':
        map = Basemap(lon_0=25,boundinglat=0,
                resolution='l',area_thresh=2000.,projection='npstere')
        levels = np.arange(-2.0, 25.0, 0.5)
        
    x, y = map(tlon,tlat)
    
    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    map.drawmapboundary()

    #map.drawmeridians(np.arange(0,360,5))
    #map.drawparallels(np.arange(0,90,5))

    temp = np.ma.masked_values(temp,grdROMS.fill_value)
    
    if var=='ssh':
        levels = np.arange(-1.5, 1.5, 0.1)
    if var in ['uvel','vvel','urot','vrot']:
        levels = np.arange(-2, 2, 0.1)
    if var=='sodamask':
        levels=np.arange(0,1,1)
    CS1 = map.contourf(x,y,temp,levels,cmap=cm.get_cmap('jet',len(levels)-1))
    #CS1 = map.contourf(x,y,temp,levels,cmap=cm.get_cmap('jet',len(levels)-1))
    plt.colorbar(orientation='horizontal')
    plt.title('Variabel %s at depth %s'%(var,depthlevel))
    #plotfile='figures/soda_depthK_'+str(depthlevel)+'.png'
    #plt.savefig(plotfile)
    plt.show()
   