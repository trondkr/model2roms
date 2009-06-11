from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile
from pylab import *

def contourMap(grdROMS,grdMODEL,data,depthlevel,var):
    
    
    tlat      = np.array(grdROMS.lat_rho)
    tlon      = np.array(grdROMS.lon_rho)
    
    if var in ['vrot','vbar','vvel']:
        tlat      = np.array(grdROMS.lat_v)
        tlon      = np.array(grdROMS.lon_v)
    if var in ['urot','ubar','uvel']:
        tlat      = np.array(grdROMS.lat_u)
        tlon      = np.array(grdROMS.lon_u)
    
    temp      = data
    
    
    plt.figure(figsize=(12,12))

    if grdROMS.grdName=='GOM':
    # Plot the Georges Bank region
        map = Basemap(llcrnrlon=-78.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
                      resolution='i',projection='tmerc',lon_0=-70,lat_0=0,area_thresh=10.)
        levels = np.arange(0.0, 26.0, 0.5)
        
    if grdROMS.grdName=='Nordic' or grdROMS.grdName=='Nordic2':
    # Plot the Nordic region (Greenland, Nordic Seas, and the Barents Sea)
        map = Basemap(lon_0=25,boundinglat=50,
                      resolution='l',area_thresh=100.,projection='npstere')
        levels = np.arange(-2.0, 15.0, 0.1)
        
    if grdROMS.grdName=='NA':
        map = Basemap(lon_0=25,boundinglat=0,
                resolution='l',area_thresh=2000.,projection='npstere')
        levels = np.arange(-2.0, 25.0, 0.5)
    
    if grdROMS.grdName=='NA_Nathan':
        map = Basemap(lon_0=25,boundinglat=0,
                resolution='l',area_thresh=2000.,projection='npstere')
        levels = np.arange(2.0, 32.0, 0.5)
        
    if grdROMS.grdName=='GOM_Nathan':
        map = Basemap(llcrnrlon=grdROMS.lon_rho[0,:].min()-0.25,
                      llcrnrlat=grdROMS.lat_rho[:,0].min()-2.25,
                      urcrnrlon=grdROMS.lon_rho[0,:].max()+5.25,
                      urcrnrlat=grdROMS.lat_rho[:,0].max()+0.25,
                      resolution='i',projection='tmerc',lon_0=-80,lat_0=0,area_thresh=10.)
        levels = np.arange(2.0, 40.0, 0.5)
        
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
    if var in ['uvel','vvel','urot','vrot','ubar','vbar']:
        levels = np.arange(-2.0, 2.0, 0.1)
    if var=='sodamask':
        levels=np.arange(0,1,1)
    CS1 = map.contourf(x,y,temp,levels,cmap=cm.get_cmap('jet',len(levels)-1) )#,alpha=0.5)
    #CS2 = contour(x,y,temp,CS1.levels[::2],
    #                    colors = 'k',
    #                    hold='on')

    plt.colorbar(orientation='horizontal')
    plt.title('Variabel %s at depth %s'%(var,depthlevel))
    #plotfile='figures/Mexico_depthK_'+str(depthlevel)+'.png'
    #plt.savefig(plotfile)
    plt.show()
   
