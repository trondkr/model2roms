
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

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2008, 8, 19)
__modified__ = datetime.datetime(2009, 1, 22)
__version__  = "1.0"
__status__   = "Development"


def getStationTime(grdSODA,year,ID):
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

def getStationData(years,IDS,outfilename,sodapath,latlist,lonlist):
      
    fileNameIn=sodapath+'SODA_2.0.2_'+str(years[0])+'_'+str(IDS[0])+'.cdf'
    
    """First time in loop, get the essential old grid information"""
    """SODA data already at Z-levels. No need to interpolate to fixed depths, but we use the one we have"""
    
    grdSODA = grd.grdClass(fileNameIn,"SODA")
    IOverticalGrid.get_z_levels(grdSODA)
    
    station=0
     
    for lat, lon in zip(latlist, lonlist):
        """Now we want to find the indices for our longitude, latitude station pairs in the lat-long list"""
        lonindex,latindex=getStationIndices(grdSODA,lon,lat)
        index=(len(years)*len(IDS),len(grdSODA.depth))
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
  
                jdsoda, yyyymmdd = getStationTime(grdSODA,year,ID)
                stTime.append(jdsoda)
                stDate.append(yyyymmdd)
                
                cdf = Dataset(filename)

                """Each SODA file consist only of one time step. Get the subset data selected, and
                store that time step in a new array:"""
                stTemp[time,:] = np.array(cdf.variables["TEMP"][0,:,latindex,lonindex])
                stSalt[time,:] = np.array(cdf.variables["SALT"][0,:,latindex,lonindex])
                stSSH[time]    = np.array(cdf.variables["SSH"][0,latindex,lonindex])
                stUvel[time,:] = np.array(cdf.variables["U"][0,:,latindex,lonindex])
                stVvel[time,:] = np.array(cdf.variables["V"][0,:,latindex,lonindex])
                
                cdf.close()
                    
                time+=1
         
        print 'Total time steps saved to file %s for station %s'%(time,station)
        plotStation.contourData(stTemp,stTime,stDate,grdSODA.depth)
        
        outfilename='Station_st%i_%s_to_%s.nc'%(station+1,years[0],years[-1])
        print 'Results saved to file %s'%(outfilename)
        writeStation(stTemp,stSalt,stUvel,stVvel,stSSH,stTime,grdSODA.depth,lat,lon,outfilename)
        station+=1
    
def getStationIndices(grdSODA,st_lon,st_lat):
    """
    This is a function that takes longitude and latitude as
    decimal input, and returns the index values closest to
    the longitude and latitude. This is an iterative process for finding the best
    index pair. Trond Kristiansen, 11.03.2008
    """
    if st_lon<0: st_lon=st_lon+360; NEG=True
    else: NEG=False
    
    longitude=grdSODA.lon
    latitude =grdSODA.lat 
    #Find closest latitude column 
    tmp2=np.abs(np.subtract(latitude,st_lat))

    #Sort the column to find the index that closest matches station
    delta_lat = tmp2.argsort(axis=0)
    
    # Find closest longitude index in latitude column
    tmp1=np.abs(np.subtract(longitude[delta_lat[0,0],:],st_lon))
    
    # Sort the result to find the best longitude in that column
    delta_lon = tmp1.argsort()
    
    tmp2=np.abs(np.subtract(latitude[:,delta_lon[0]],st_lat))

    # Repeat process to find the best latitude in the best longitude column
    delta_lat = tmp2.argsort(axis=0)
    
    tmp1=np.abs(np.subtract(longitude[delta_lat[0],:],st_lon))
    
    delta_lon = tmp1.argsort()
    
   
    print '' 
    print '=====get_station======'
    if NEG is True:
        print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon-360,st_lat)
    else:
        print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon,st_lat)
    print 'Result ===>'
    print 'Found index pair [%s,%s] in gridfile'%(delta_lon[0],delta_lat[0])
    if NEG is True:
        print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[delta_lat[0],delta_lon[0]]-360.,latitude[delta_lat[0],delta_lon[0]])
    else:
        print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[delta_lat[0],delta_lon[0]],latitude[delta_lat[0],delta_lon[0]])
    print '======================'
    print ''
        
    return delta_lon[0], delta_lat[0]



def writeStation(t,s,uvel,vvel,ssh,ntime,depth,lat,lon,outfilename):
       
    if os.path.exists(outfilename):
        os.remove(outfilename)
    print 'Writing first results to file %s'%(outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF3_CLASSIC')
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a station (lat=%3.2f,lon=%3.2f) file from SODA data"%(lat,lon)
    f1.history = 'Created ' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 time-depth file created using SODA2ROMS' 
       
    # Define dimensions
    f1.createDimension('time', len(ntime))
    f1.createDimension('z', len(depth))
     
    v=f1.createVariable('lon','d',zlib=True)
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
 
    
    f1.close()



