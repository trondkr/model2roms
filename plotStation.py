from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, NetCDFFile
from pylab import *

def contourData(data,timedata,datedata,depthdata):
    
    plt.figure(figsize=(12,12))

    x = timedata
    y = -flipud(depthdata)

    X,Y = meshgrid(x, y)
    print datedata
     
    pcolor(X,Y,rot90(data),shading = "flat")
    plt.colorbar(orientation='horizontal')
   
    axis('tight')
    clim([0,25]) 
    ylim([-500,-5])
    plt.title('test of contour')
    
    sep=73
    s=0
    d=[]
    t=[]
    
    for i in range(len(timedata)):
        if s==sep:
            d.append(datedata[i])
            t.append(timedata[i])
            s=0
           
        s+=1
    locs, labels = xticks(t, d)
    setp(labels, 'rotation', 'vertical')
    
    show()

 
    
    #plotfile='figures/soda_depthK_'+str(depthlevel)+'.png'
    #plt.savefig(plotfile)