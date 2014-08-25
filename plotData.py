import os
from mpl_toolkits.basemap import Basemap
from pylab import *

def contourMap(grdROMS, tlon, tlat, mydata, depthlevel, var):
    print "NAME:",grdROMS.grdName
    map = Basemap(projection='ortho',lon_0=-20,lat_0=25, resolution='c', area_thresh=10000)
    levels = np.arange(-2.0, 30.0, 0.5)

    plt.figure(figsize=(10,10), frameon=False)
    #mydata = np.ma.masked_where(mydata==grdROMS.fill_value,mydata)

    if grdROMS.grdName=='GOM':
    # Plot the Georges Bank region
        map = Basemap(llcrnrlon=-78.5,llcrnrlat=32.5,urcrnrlon=-58,urcrnrlat=45.5,
                      resolution='i',projection='tmerc',lon_0=-70,lat_0=0,area_thresh=10.)
        levels = np.arange(0.0, 26.0, 0.5)

    if grdROMS.grdName=='Nordic' or grdROMS.grdName=='Nordic2':
        map = Basemap(lon_0=25,boundinglat=50,
                      resolution='l',area_thresh=100.,projection='npstere')
        levels = np.arange(-2.0, 20.0, 0.1)

    if grdROMS.grdName=='Troms':
        map = Basemap(projection='merc',llcrnrlat=68,urcrnrlat=72,\
            llcrnrlon=15,urcrnrlon=24,lat_ts=70,resolution='h',area_thresh=0.)
        levels = np.arange(-2.0, 8.0, 0.25)

    if grdROMS.grdName=='NS8KM':
        map = Basemap(llcrnrlon=-18.5,
                      llcrnrlat=43.5,
                      urcrnrlon=32.5,
                      urcrnrlat=67.7,
                      resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=20.)
        levels = np.arange(-2.0, 20.0, 0.5)

    if grdROMS.grdName=='NA':
        map = Basemap(lon_0=25,boundinglat=0,
                resolution='l',area_thresh=2000.,projection='npstere')
        levels = np.arange(-2.0, 25.0, 0.5)

    if grdROMS.grdName=='NATL':
        map = Basemap(lon_0=2,boundinglat=0,
                resolution='l',area_thresh=2000.,projection='npstere')
        levels = np.arange(-2.0, 30.0, 0.5)

    if grdROMS.grdName=='REGSCEN':
        print "Plotting REGSCEN"
        map = Basemap(projection='ortho',lon_0=-20,lat_0=25, resolution='c', area_thresh=10000)
        levels = np.arange(-2.0, 30.0, 0.5)
    print "Test",grdROMS.grdName
    x, y = map(tlon,tlat)

    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()

    if var=='temperature':
        levels = np.arange(-2.0, 36.0, 1)
    if var=='salinity':
        levels = np.arange(15.0, 40.0, 0.3)
    if var=='runoff':
        levels = np.arange(1e-4, 1e-8, 1e-10)
    if var=='ssh':
        levels = np.arange(-1.5, 1.5, 0.1)
    if var in ['uvel','vvel','urot','vrot','ubar','vbar']:
        levels = np.arange(-1.0, 1.0, 0.1)
    if var=='sodamask':
        levels=np.arange(0,1,1)
    if var=='hice':
        levels=np.arange(0,100,0.1)
    if var=='aice':
        levels=np.arange(0,100,5)
    if var=='uice' or var=="vice":
        levels=np.arange(-100,100,5)
    if var=='ageice':
        levels=np.arange(0,12,1)
    if var=='snow_thick':
        levels=np.arange(0,2,0.1)
    if var =='runoff':
        CS1 = map.contourf(x, y, mydata)#,alpha=0.5)
        plt.colorbar(orientation='horizontal')
    else:
        CS1 = map.contourf(x, y, mydata, levels, cmap=cm.get_cmap('RdYlBu_r',len(levels)-1) )#,alpha=0.5)
        plt.colorbar(CS1, orientation='vertical', extend='both', shrink=0.5)

   
    plt.title('Var:%s - depth:%s - time:%s'%(var,depthlevel,grdROMS.time))
    plotfile='figures/'+str(var)+'_depth_ESMF_'+str(depthlevel)+'_time_'+str(grdROMS.time)+'.png'
    if not os.path.exists('figures'):
        os.makedirs('figure')
    plt.savefig(plotfile)
    print "Saved figure: %s"%(plotfile)
    #plt.show()

def contourStationData(data,timedata,datedata,depthdata,stationName):

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
