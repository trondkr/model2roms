
import os, sys, datetime
import numpy as np
import grd
import IOverticalGrid
import soda2roms
import date
from netCDF4 import Dataset
from netCDF4 import num2date
import time
import os
import plotStation
import Nio

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2008, 8, 19)
__modified__ = datetime.datetime(2009, 1, 22)
__version__  = "1.0"
__status__   = "Development"

def getAverage(latindex,lonindex,index):
    """The input file for this routine is created by making an average of all the SODA files from 1961-1990.
    This is done using the script soda2average/average_soda.py"""
    
    averageFile='../soda2average/average_SODA_1961-1990.nc'
    if os.path.exists(averageFile):
        print 'Average file %s exists so anomaly values of salt, '%(averageFile)
        print 'temp, u, and v will be written to file'
        ave = Dataset(averageFile,'r')
        
        average=True
        aveTemp=np.zeros((index),np.float64)
        aveSalt=np.zeros((index),np.float64)
        aveUvel=np.zeros((index),np.float64)
        aveVvel=np.zeros((index),np.float64)
    
        aveTemp[:] = ave.variables["temp"][:,latindex,lonindex]
        aveSalt[:] = ave.variables["salt"][:,latindex,lonindex]
        aveUvel[:] = ave.variables["uvel"][:,latindex,lonindex]
        aveVvel[:] = ave.variables["vvel"][:,latindex,lonindex]
        ave.close()
    else:
        average=False
        aveTemp=None; aveSalt=None; aveUvel=None; aveVvel=None
        
    return aveTemp,aveSalt,aveUvel,aveVvel, average

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
    
    print 'Current time of SODA file : %s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    
    return jdsoda-jdref, yyyymmdd

def getStationData(years,IDS,sodapath,latlist,lonlist):
      
    fileNameIn=sodapath+'SODA_2.0.2_'+str(years[0])+'_'+str(IDS[0])+'.cdf'
    
    """First time in loop, get the essential old grid information"""
    """SODA data already at Z-levels. No need to interpolate to fixed depths, but we use the one we have"""
    
    grdMODEL = grd.grdClass(fileNameIn,"SODA")
    IOverticalGrid.get_z_levels(grdMODEL)
    
    station=0
    numberOfPoints=4
    
    for lat, lon in zip(latlist, lonlist):
        """Now we want to find the indices for our longitude, latitude station pairs in the lat-long list"""
        gridIndexes, dis = getStationIndices(grdMODEL,lon,lat,'SODA',numberOfPoints)
        
        index=(len(years)*len(IDS),len(grdMODEL.depth))
        indexSSH=(len(years)*len(IDS))
        
        print '-------------------------------------------------------------\n'
        print 'Extracting data for station (%s,%s) for the years %s->%s'%(lon,lat,years[0],years[-1])
        print 'Size of output data array is ->',index
        
        stTemp=np.zeros((index),np.float64)
        stSalt=np.zeros((index),np.float64)
        stSSH =np.zeros((indexSSH),np.float64)
        stUvel=np.zeros((index),np.float64)
        stVvel=np.zeros((index),np.float64)
        stTime=[]
        stDate=[]
        time=0
       
        for year in years:
            for ID in IDS:
                file="SODA_2.0.2_"+str(year)+"_"+str(ID)+".cdf"
                filename=sodapath+file
               
                jdsoda, yyyymmdd = getStationTime(grdMODEL,year,ID)
                stTime.append(jdsoda)
                stDate.append(yyyymmdd)
                t=0
                cdf = Dataset(filename,'r',format='NETCDF3')
 
                """Each SODA file consist only of one time step. Get the subset data selected, and
                store that time step in a new array:"""
                for i in range(numberOfPoints):
                    wgt=float(dis[i])/sum(dis)
                    latindex=int(gridIndexes[i][0])
                    lonindex=int(gridIndexes[i][1])

                    """The values at a station is calculated by interpolating from the
                    numberOfPoints around the station uysing weights (wgt)
                    """
                    stTemp[time,:] = stTemp[time,:] + (cdf.variables["TEMP"][t,:,latindex,lonindex])*wgt
                    stSalt[time,:] = stSalt[time,:] + (cdf.variables["SALT"][t,:,latindex,lonindex])*wgt
                    stSSH[time]    = stSSH[time]    + (cdf.variables["SSH"][t,latindex,lonindex])*wgt
                    stUvel[time,:] = stUvel[time,:] + (cdf.variables["U"][t,:,latindex,lonindex])*wgt
                    stVvel[time,:] = stVvel[time,:] + (cdf.variables["V"][t,:,latindex,lonindex])*wgt
                 
                cdf.close()
                    
                time+=1
         
        print 'Total time steps saved to file %s for station %s'%(time,station)
        #plotStation.contourData(stTemp,stTime,stDate,grdMODEL.depth)
        
        aveTemp,aveSalt,aveUvel,aveVvel,average = getAverage(latindex,lonindex,len(grdMODEL.depth))
        
        outfilename='Station_st%i_%s_to_%s.nc'%(station+1,years[0],years[-1])
        outfilename='station_lat_'+str(lat)+'_lon_'+str(lon)+'.xyz'
        print 'Results saved to file %s'%(outfilename)
        writeStationNETCDF4(stTemp,stSalt,stUvel,stVvel,stSSH,stTime,grdMODEL.depth,lat,lon,outfilename,average,aveTemp,aveSalt,aveUvel,aveVvel)
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
    
    if type=='SODA':
        longitude=grdObject.lon
        latitude =grdObject.lat
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
    print '=====get_station======'
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



def writeStationNETCDF4(t,s,uvel,vvel,ssh,ntime,depth,lat,lon,outfilename,average,aveTemp,aveSalt,aveUvel,aveVvel):
       
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
    
    v_u=f1.createVariable('u', 'f', ('time','z'), zlib=True)
    v_u.long_name = "U-velocity, scalar, series"
    v_u.units = "m/s"
    v_u[:,:] = uvel
    
    v_v=f1.createVariable('v', 'f', ('time','z'), zlib=True)
    v_v.long_name = "V-velocity, scalar, series"
    v_v.units = "m/s"
    v_v[:,:] = vvel
 
    if average==True:
        v_anom=f1.createVariable('tempAnomaly', 'f', ('time','z'), zlib=True)
        v_anom.long_name = "Ocean temperature anomaly"
        v_anom.units = "degrees Celsius"

        for k in range(len(t[:,1])):
             v_anom[k,:] = t[k,:] - aveTemp[:]
       
        v_anom=f1.createVariable('saltAnomaly', 'f', ('time','z'), zlib=True)
        v_anom.long_name = "Ocean salinity anomaly"
        v_anom.units = "psu"

        for k in range(len(s[:,1])):
             v_anom[k,:] = s[k,:] - aveSalt[:]
    
        v_anom=f1.createVariable('uAnomaly', 'f', ('time','z'), zlib=True)
        v_anom.long_name = "Ocean u anomaly"
        v_anom.units = "m/s"

        for k in range(len(s[:,1])):
             v_anom[k,:] = uvel[k,:]- aveUvel[:]
     
     
        v_anom=f1.createVariable('vAnomaly', 'f', ('time','z'), zlib=True)
        v_anom.long_name = "Ocean v anomaly"
        v_anom.units = "m/s"

        for k in range(len(s[:,1])):
             v_anom[k,:] = vvel[k,:]- aveVvel[:]
             
    f1.close()


def writeStationNETCDF3(t,s,uvel,vvel,ssh,ntime,depth,lat,lon,outfilename):
       
    if os.path.exists(outfilename):
        os.remove(outfilename)
    print 'Writing first results to file %s'%(outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
    f1.description="This is a station (lat=%3.2f,lon=%3.2f) file from SODA data"%(lat,lon)
    f1.history = 'Created ' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF3 time-depth file created using SODA2ROMS' 
       
    # Define dimensions
    f1.createDimension('time', len(ntime))
    f1.createDimension('z', len(depth))
     
    v=f1.createVariable('lon','d')
    v.long_name = "Longitude position of station" ;
    v.units = "Longitude"
    v[:]=lon
    
    v=f1.createVariable('lat','d')
    v.long_name = "Latitude position of station" ;
    v.units = "Latitude"
    v[:]=lat
    
    v=f1.createVariable('depth','d',('z',))
    v.long_name = "Z-depth matrix" ;
    v.units = "meter"
    v[:]=depth
     
    v_time = f1.createVariable('time', 'd', ('time',))
    v_time.long_name = 'Days since 1948-01-01 00:00:00'
    v_time.units = 'days'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = ntime
    
    v_temp=f1.createVariable('temp', 'f', ('time','z'))
    v_temp.long_name = "Ocean temperature"
    v_temp.units = "degrees Celsius"
    v_temp[:,:] = t
    
    v_salt=f1.createVariable('salt', 'f', ('time','z'))
    v_salt.long_name = "Ocean salinity"
    v_salt.units = "psu"
    v_salt[:,:]=s
    
    v_ssh=f1.createVariable('zeta','d',('time',))
    v_ssh.long_name = "Sea surface height (SSH)"
    v_ssh.units = "m"
    v_ssh[:] = ssh
    
    v_u=f1.createVariable('u', 'f', ('time','z'))
    v_u.long_name = "U-velocity, scalar, series"
    v_u.units = "m/s"
    v_u[:,:] = uvel
    
    v_v=f1.createVariable('v', 'f', ('time','z'))
    v_v.long_name = "V-velocity, scalar, series"
    v_v.units = "m/s"
    v_v[:,:] = vvel
 
    
    f1.close()

