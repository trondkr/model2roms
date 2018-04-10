import os
from mpl_toolkits.basemap import Basemap
from pylab import *

def contourMap(grdROMS, tlon, tlat, mydata1, mydata2, mydata3, var, mytype, currentdate):

    plt.figure(figsize=(10,10), frameon=False)
    
    map = Basemap(lon_0=25,boundinglat=50,
                      resolution='l',area_thresh=100.,projection='npstere')
    
    x, y = list(map(tlon,tlat))

    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()

    if var=='wind':
        levels = np.arange(np.min(mydata3),np.max(mydata3),0.1)
   
        CS1 = map.contourf(x, y, mydata3, levels, cmap=cm.get_cmap('RdYlBu_r',len(levels)-1) )#,alpha=0.5)
        plt.colorbar(CS1, orientation='vertical', extend='both', shrink=0.5)
            
        if mytype=="REGSCEN":
          step=8
        else:
          step=1
        map.quiver(x[0:-1:step,0:-1:step],y[0:-1:step,0:-1:step],
            mydata1[0:-1:step,0:-1:step],mydata2[0:-1:step,0:-1:step],
            scale=400)
   
  #  plt.title('Var:%s - depth:%s - time:%s'%(var,grdROMS.time))
    plotfile='figures/'+str(var)+'_'+str(mytype)+'_time_'+str(currentdate)+'.png'
    if not os.path.exists('figures'):
        os.makedirs('figure')
    plt.savefig(plotfile)
    print("Saved figure: %s"%(plotfile))
   # plt.show()
