from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile

def contourMap(grd,data):
    
 
    tlat      = np.array(grd.lat_rho)
    tlon      = np.array(grd.lon_rho)
    # masked array returned, masked where data == _FillValue
    temp      = data
    
    
    # make longitudes monotonically increasing.
    tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
    
    plt.figure(figsize=(6,8))
    plt.subplot(2,1,1)
    # Plot the Georges Bank region
    #map = Basemap(llcrnrlon=-80.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
    #        resolution='i',projection='tmerc',lon_0=-70,lat_0=0)
    
    # Plot the Nordic region
    map = Basemap(llcrnrlon=-20,llcrnrlat=40.5,urcrnrlon=70,urcrnrlat=80,
            resolution='c',projection='tmerc',lon_0=10,lat_0=60)
    
    map.drawcoastlines()
    map.fillcontinents(color='white')
    
    x, y = map(tlon,tlat)
    im = map.pcolor(x,y,temp,\
               shading='faceted',cmap=plt.cm.cool,vmin=0,vmax=0)
    # disclaimer:  these are not really the grid cells because of the
    # way pcolor interprets the x and y args.
    plt.title('(A) ROMSGrid Cells')
    
    # subplot 2 is a contour plot of surface temperature from the
    # CCSM ocean model.
    plt.subplot(2,1,2)
    map.drawcoastlines()
    map.fillcontinents(color='white')
    
    CS1 = map.contourf(x,y,temp,15)
    #CS2 = map.contour(x,y,temp,15,colors='black',linewidths=0.5)
    plt.title('(B) Surface Temp contours on POP Grid')
    
    plt.show()
    #plt.savefig('roms_grid.png')