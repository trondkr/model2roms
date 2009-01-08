from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile

def contourMap(grdROMS,grdSODA,data,depthlevel,var):
    
 
    tlat      = np.array(grdROMS.lat_rho)
    tlon      = np.array(grdROMS.lon_rho)
    temp      = data
    
    
    # make longitudes monotonically increasing.
    tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
    
    plt.figure(figsize=(12,12))

    # Plot the Georges Bank region
    #map = Basemap(llcrnrlon=-80.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
    #        resolution='i',projection='tmerc',lon_0=-70,lat_0=0)
    
    # Plot the Nordic region
    map = Basemap(llcrnrlon=-80,llcrnrlat=31.5,urcrnrlon=-55,urcrnrlat=45,
            resolution='h',projection='tmerc',lon_0=-65,lat_0=40)
   
    map.drawcoastlines()
    map.fillcontinents(color='white')
    
    x, y = map(tlon,tlat)
    
    x2,y2=map(np.array(grdSODA.lon),np.array(grdSODA.lat))
   # im = map.pcolor(x,y,temp,\
   #            shading='faceted',cmap=plt.cm.cool,vmin=0,vmax=0)
    # disclaimer:  these are not really the grid cells because of the
    # way pcolor interprets the x and y args.
    #plt.title('Variabel % at depth %s'%(var,depthlevel))
    
    # subplot 2 is a contour plot of surface temperature from the
    # CCSM ocean model.
  
    map.drawcoastlines()
    map.fillcontinents(color='white')
    
    CS1 = map.contourf(x,y,temp,15)
    plt.colorbar(orientation='horizontal')
    #CS2 = map.contour(x,y,temp,15,colors='black',linewidths=0.5)
    plt.title('Variabel %s at depth %s'%(var,depthlevel))
    
    map.plot(x2,y2,'bo')
     
    plt.show()
    #plotfile='figures/soda_depthK_'+str(depthlevel)+'.png'
    #plt.savefig(plotfile)