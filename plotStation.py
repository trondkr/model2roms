from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile
from pylab import *

def contourData(data,timedata,datedata,depthdata,stationName):
    
    x = timedata
    y = flipud(depthdata)
    X,Y = meshgrid(x, y)
    data=rot90(data)
    levels=np.arange(0, data.max(), 0.5)
    plt.figure(figsize=(10,4))
    cs1=contourf(X,Y,data,levels,cmap=cm.get_cmap('RdYlBu_r',len(levels)-1),alpha=1.0,extend='both')
    cs2=contour(X,Y,data,cs1.levels[::4],colors = ('k',),hold='on',linewidths = (0.1,))
    colorbar(cs1,extend='both')
   
    sep=73; s=0; d=[]; t=[]
    
    for i in range(len(timedata)):
        if s==sep:
            d.append(datedata[i])
            t.append(timedata[i])
            s=0
           
        s+=1
    locs, labels = plt.xticks(t, d)
    plt.setp(labels, 'rotation', 'vertical')
    
    plotfile='figures/station_'+str(stationName)+'.png'
    plt.savefig(plotfile)
    
    #plt.show()