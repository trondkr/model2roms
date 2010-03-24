
import os, sys, datetime
import numpy as np
import grd
import IOverticalGrid
import soda2roms
import date
from netCDF4 import Dataset
from netCDF4 import num2date
import os, time
import plotData
import string
from progressBar import progressBar

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2008, 8, 19)
__modified__ = datetime.datetime(2009, 1, 22)
__version__  = "1.0"
__status__   = "Development"

def getAverage(yyyymmdd,gridIndexes,validIndex,validDis):
    """The input file for this routine is created by making an average of all the SODA files from 1961-1990.
    This is done using the script soda2average/averageSODA.py"""
    """First get month of interest:"""
    d = string.split(yyyymmdd,'/'); month=int(d[1])
    averageFile='../soda2average/clim/averageSODA1961-1990.nc'
    if os.path.exists(averageFile):
        print '======================================================================='
        print 'Average file %s exists so clim values written to file'%(averageFile)  
        print '======================================================================='
        ave = Dataset(averageFile,'r')
        
        index =(len(ave.variables["depth"][:]))
        average=True
        aveTemp=np.zeros((index),np.float64)
        aveSalt=np.zeros((index),np.float64)
        aveUvel=np.zeros((index),np.float64)
        aveVvel=np.zeros((index),np.float64)
        
        for i in validIndex:
            wgt=float(validDis[i])/sum(validDis)
            latindex=int(gridIndexes[i][0])
            lonindex=int(gridIndexes[i][1])
            
            """The values at a station is calculated by interpolating from the
            numberOfPoints around the station uysing weights (wgt)
            """
            print "test",i, month, latindex,lonindex, ave.variables["temp"].shape
            aveTemp[:] = aveTemp[:]  + (ave.variables["temp"][:,latindex,lonindex,month-1])*wgt
            aveSalt[:] = aveSalt[:]  + (ave.variables["salt"][:,latindex,lonindex,month-1])*wgt
            aveUvel[:] = aveUvel[:]  + (ave.variables["uvel"][:,latindex,lonindex,month-1])*wgt
            aveVvel[:] = aveVvel[:]  + (ave.variables["vvel"][:,latindex,lonindex,month-1])*wgt
         
        ave.close()
        print aveTemp, aveSalt, aveUvel, aveVvel
    else:
        average=False
        aveTemp=None; aveSalt=None; aveUvel=None; aveVvel=None; aveTime=None
        
    return aveTemp,aveSalt,aveUvel,aveVvel,average

def getStationTime(grdMODEL,year,ID):
    """
    Find the day and month that the SODA file respresents based on the year and ID number.
    Each SODA file represents a 5 day average, therefore we let the date we find be the first day
    of those 5 days. Thats the reason we subtract 4 below for day of month.
    """
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    """
    ref_date = date.Date()
    ref_date.day=1
    ref_date.month=1
    ref_date.year=1948
    jdref=ref_date.ToJDNumber()
    
    days=0.0; month=1;loop=True
    
    while loop is True:
        
        d=date.NumberDaysMonth(month,year)
        if days+d<int(ID)*5:
            days=days+d
            month+=1
        else:
            day=int(int(ID)*5-days)
            loop=False
            
    soda_date = date.Date()
    soda_date.day=day
    soda_date.month=month
    soda_date.year=year
    jdsoda=soda_date.ToJDNumber()
    yyyymmdd = '%s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    
    message='Current time of SODA file : %s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    
    return jdsoda-jdref, yyyymmdd, message

def testValidStation(cdf,dis,numberOfPoints,gridIndexes):
    validIndex=[]; validDis=[]
    for i in range(numberOfPoints):
        latindex=int(gridIndexes[i][0])
        lonindex=int(gridIndexes[i][1])
        
        """Test to see if station contains anything but missing values. If it does, we assume
        this holds for salinity, u, and v as well. We also assume all values less than 10000"""
        if any(abs(cdf.variables["TEMP"][0,:,latindex,lonindex])< 10000):
            validIndex.append(i)
            validDis.append(dis[i])
  
    if not validIndex:
        print 'No valid data found for position'
        exit()
    print "Found %s valid surounding grid cells for station"%(len(validIndex))
    return validIndex, validDis

def testValidDepth(cdf,numberOfPoints,gridIndexes,depth):
    """To avoid having to extract all depth layers that do not
    contain valid data, we find the deepest depth layer of the surrouding stations
    and use that indices in the further extration"""
    deepestMax=40 # The deepes soda depth is 40
    deepest=0
    for i in range(numberOfPoints):
        latindex=int(gridIndexes[i][0])
        lonindex=int(gridIndexes[i][1])
        
        """Test to see how many depth layers we should extract. We do not want layers
        that only contain missing values"""
        d = cdf.variables["TEMP"][0,:,latindex,lonindex]
        for k in range(len(d)):
            if abs(d[k]) < 1000:
                """One of the other stations may have values valid at
                deeper stations than this one. Make sure we find absolute deepest"""
                #if k>deepest:
                deepest=k
       
    print "Found deepest valid depth-layer %s which is equivavelnt to %sm\n" %(deepest+1, depth[deepest+1])
    """Now we return deepest + 1 as the index is one value more than the counter"""
    return deepest      
    
def initArrays(years,IDS,deepest,name,lon,lat):
    index=(len(years)*len(IDS),deepest)
    indexSSH=(len(years)*len(IDS))
        
    print 'Extracting data for station (%s,%s) for the years %s->%s'%(lon,lat,years[0],years[-1])
    print 'Size of output data array is ->',index
   
    stTemp=np.zeros((index),np.float64)
    stSalt=np.zeros((index),np.float64)
    stSSH =np.zeros((indexSSH),np.float64)
    stUvel=np.zeros((index),np.float64)
    stVvel=np.zeros((index),np.float64)
    stTauX=np.zeros((indexSSH),np.float64)
    stTauY=np.zeros((indexSSH),np.float64)
    
    return stTemp, stSalt, stSSH, stUvel, stVvel, stTauX, stTauY

def getStationData(years,IDS,sodapath,latlist,lonlist,stationNames):
    empty  =u'\u25FD'; filled =u'\u25FE'
    progress  = progressBar(color='red',width=30, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))
    
    fileNameIn=sodapath+'SODA_2.0.2_'+str(years[0])+'_'+str(IDS[0])+'.cdf'
    
    """First time in loop, get the essential old grid information"""
    """SODA data already at Z-levels. No need to interpolate to fixed depths, but we use the one we have"""
    
    grdMODEL = grd.grdClass(fileNameIn,"SODA")
    IOverticalGrid.get_z_levels(grdMODEL)
    
    station=0
    numberOfPoints=4
   
    for lat, lon in zip(latlist, lonlist):
        print '\n----------------NEW STATION==> %s ------------------------------------------'%(stationNames[station])
     
        """Now we want to find the indices for our longitude, latitude station pairs in the lat-long list"""
        gridIndexes, dis = getStationIndices(grdMODEL,lon,lat,'SODA',numberOfPoints)
    
        stTime=[]; stDate=[]; time=0; counter=0; t=0
        total=float(len(years)*len(IDS))
        for year in years:
            for ID in IDS:
                file="SODA_2.0.2_"+str(year)+"_"+str(ID)+".cdf"
                filename=sodapath+file
               
                jdsoda, yyyymmdd, message = getStationTime(grdMODEL,year,ID)
                
                stTime.append(jdsoda)
                stDate.append(yyyymmdd)
               
                cdf = Dataset(filename,'r')
 
                """Each SODA file consist only of one time step. Get the subset data selected, and
                store that time step in a new array:"""
                if year==years[0] and ID==IDS[0]:
                    validIndex, validDis = testValidStation(cdf,dis,numberOfPoints, gridIndexes)
                    deepest              = testValidDepth(cdf,numberOfPoints, gridIndexes,grdMODEL.depth)
                    stTemp, stSalt, stSSH, stUvel, stVvel, stTauX, stTauY = initArrays(years,IDS,deepest,stationNames[station],lon,lat)
                
                for i in validIndex:
                    wgt=float(validDis[i])/sum(validDis)
                    latindex=int(gridIndexes[i][0])
                    lonindex=int(gridIndexes[i][1])
                    
                    """The values at a station is calculated by interpolating from the
                    numberOfPoints around the station uysing weights (wgt)
                    """
                    stTemp[time,:] = stTemp[time,:] + (cdf.variables["TEMP"][t,0:deepest,latindex,lonindex])*wgt
                    stSalt[time,:] = stSalt[time,:] + (cdf.variables["SALT"][t,0:deepest,latindex,lonindex])*wgt
                    stSSH[time]    = stSSH[time]    + (cdf.variables["SSH"][t,latindex,lonindex])*wgt
                    stTauX[time]    = stTauX[time]    + (cdf.variables["TAUX"][t,latindex,lonindex])*wgt
                    stTauY[time]    = stTauY[time]    + (cdf.variables["TAUY"][t,latindex,lonindex])*wgt
                    stUvel[time,:] = stUvel[time,:] + (cdf.variables["U"][t,0:deepest,latindex,lonindex])*wgt
                    stVvel[time,:] = stVvel[time,:] + (cdf.variables["V"][t,0:deepest,latindex,lonindex])*wgt
                 
                cdf.close()
                counter+=1
                
                """Show progress bar"""
                p=int( ((time*1.0+1)/(1.0*total))*100.)
                progress.render(p,message)
                time+=1
                
        
        print 'Total time steps saved to file %s for station %s'%(time,station)
        #plotData.contourStationData(stTemp,stTime,stDate,-grdMODEL.depth[0:deepest],stationNames[station])
        
        outfilename='station_'+str(stationNames[station])+'.nc'
        print 'Results saved to file %s'%(outfilename)
        writeStationNETCDF4(stTemp,stSalt,stUvel,stVvel,stSSH,stTauX,stTauY,stTime,
                            grdMODEL.depth[0:deepest],lat,lon,outfilename)
        station+=1
    
def getStationIndices(grdObject,st_lon,st_lat,type,numberOfPoints):
    """
    This is a function that takes longitude and latitude as
    decimal input, and returns the index values closest to
    the longitude and latitude. This is an iterative process for finding the best
    index pair. Trond Kristiansen, 11.03.2008, 09.06.2009
    """
    if st_lon<0: st_lon=st_lon+360.0; NEG=True
    else: NEG=False
    
    if type=='SODA' or type=='AVERAGE':
        longitude=grdObject.lon
        latitude =grdObject.lat
        """Input longitude should go from 0-360"""
        longitude=np.where(longitude<0,longitude+360,longitude)
    if type=='ROMS':
        longitude=grdObject.lon_rho
        latitude =grdObject.lat_rho 
    
    distance = np.zeros((longitude.shape),dtype=np.float64)
    listd=[]
    """First, create a list of distances from the station of interest, while
    also save the matrix of dostances that contains the info to get the index pair that the distance
    of interest corresponds to"""
    for eta in range(len(latitude[:,0])):
        for xi in range(len(latitude[0,:])):
            distance[eta,xi] = np.sqrt( (latitude[eta,xi]-st_lat)**2.0 + (longitude[eta, xi] - st_lon)**2.0 )
            listd.append(distance[eta,xi])

    listsIndexes=[]
    listd.sort()
    
    """Now find the closest point to the station. When that point is found, remove the
    closests pooint and find the next closests point, until you have found numberOfPoints
    closests to station.
    """
    for i in range(numberOfPoints):
        value=listd[0]
        itemindex=np.where(distance==value)
        listsIndexes.append(itemindex)
       
        listd.pop(0)
  
    print '' 
    print '=====getStationIndices======'
    if NEG is True:
        print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon-360,st_lat)
    else:
        print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon,st_lat)
    print 'Result ===>'
    for i in range(numberOfPoints):
        print 'Found index pair in gridfile',listsIndexes[i]
        if NEG is True:
            print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[listsIndexes[i][0],listsIndexes[i][1]]-360,latitude[listsIndexes[i][0],listsIndexes[i][1]])
        else:
            print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[listsIndexes[i][0],listsIndexes[i][1]],latitude[listsIndexes[i][0],listsIndexes[i][1]])
    print '======================'
    print ''
      
    """
    We want to use data interpolated from the 4 surrounding points to get appropriate values at station point.
    We do this by using relative weights determined by relative distance to total distance from all 4 points.
    Trond Kristiansen, 09.06.2009
    """
    dis=[]
    
    for i in range(numberOfPoints):
        dis.append(np.sqrt( (latitude[listsIndexes[i][0],listsIndexes[i][1]]-st_lat)**2.0 + (longitude[listsIndexes[i][0],listsIndexes[i][1]] - st_lon)**2.0 ))
        
    return listsIndexes, dis



def writeStationNETCDF4(t,s,uvel,vvel,ssh,taux,tauy,ntime,depth,lat,lon,outfilename):
       
    if os.path.exists(outfilename):
        os.remove(outfilename)
    print 'Writing first results to file %s'%(outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a station (lat=%3.2f,lon=%3.2f) file from SODA data"%(lat,lon)
    f1.history = 'Created ' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 time-depth file created using SODA2ROMS' 
       
    f1.createDimension('time', len(ntime))
    f1.createDimension('z', len(depth))
  
    v=f1.createVariable('lon','d')
    v.long_name = "Longitude position of station" ;
    v.units = "Longitude"
    v[:]=lon
    
    v=f1.createVariable('lat','d', zlib=True)
    v.long_name = "Latitude position of station" ;
    v.units = "Latitude"
    v[:]=lat
    
    v=f1.createVariable('depth','d',('z',), zlib=True)
    v.long_name = "Z-depth matrix" ;
    v.units = "meter"
    v[:]=depth
     
    v_time = f1.createVariable('time', 'd', ('time',), zlib=True)
    v_time.long_name = 'Days since 1948-01-01 00:00:00'
    v_time.units = 'days'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = ntime
    
    v_temp=f1.createVariable('temp', 'f', ('time','z'), zlib=True)
    v_temp.long_name = "Ocean temperature"
    v_temp.units = "degrees Celsius"
    v_temp[:,:] = t
    
    v_salt=f1.createVariable('salt', 'f', ('time','z'), zlib=True)
    v_salt.long_name = "Ocean salinity"
    v_salt.units = "psu"
    v_salt[:,:]=s
    
    v_ssh=f1.createVariable('zeta','d',('time',), zlib=True)
    v_ssh.long_name = "Sea surface height (SSH)"
    v_ssh.units = "m"
    v_ssh[:] = ssh
    
    v_taux=f1.createVariable('taux','d',('time',), zlib=True)
    v_taux.long_name = "Wind stress in x direction"
    v_taux.units = "N/m^2"
    v_taux[:] = taux
    
    v_tauy=f1.createVariable('tauy','d',('time',), zlib=True)
    v_tauy.long_name = "Wind stress in y direction"
    v_tauy.units = "N/m^2"
    v_tauy[:] = tauy
    
    v_u=f1.createVariable('u', 'f', ('time','z'), zlib=True)
    v_u.long_name = "U-velocity, scalar, series"
    v_u.units = "m/s"
    v_u[:,:] = uvel
    
    v_v=f1.createVariable('v', 'f', ('time','z'), zlib=True)
    v_v.long_name = "V-velocity, scalar, series"
    v_v.units = "m/s"
    v_v[:,:] = vvel
        
    f1.close()