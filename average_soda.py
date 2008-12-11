import os, sys, datetime
import Nio
import numpy
import IOnetcdf
import plot_soda


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2008, 8, 19)
__version__  = "1.0"
__status__   = "Development"

def write_netcdf(lon, lat, data, ncFileName):
    ncFileName=ncFileName+'.cdf'
    if os.path.exists(ncFileName): os.remove(ncFileName)
    outNCFile = Nio.open_file(ncFileName, "c") 
    outNCFile.create_dimension("lon", len(lon)) 
    outNCFile.create_dimension("lat", len(lat)) 
    
    outNCFile.create_variable("lon", 'd', ('lon',)) 
    outNCFile.variables['lon'].units = "degrees_east" 
    outNCFile.variables['lon'].long_name = "Longitude" 
    outNCFile.create_variable('lat', 'd', ('lat',)) 
    outNCFile.variables['lat'].units = "degrees_north" 
    outNCFile.variables['lat'].long_name = "Latitude" 
    #outNCFile.create_variable('level', 'f', ('level',)) 
    #outNCFile.variables['level'].units = "millibar" 
    #outNCFile.variables['level'].positive = "down" 
    #outNCFile.variables['level'].long_name = "Level" 
    #outNCFile.variables['level'].desc = "Level here actually describes an ensemble forescst identifer, (see levids variable)" 
    #outNCFile.create_variable('time', 'd', ('time',)) 
    #outNCFile.variables['time'].units = "hours since 1-1-1 00:00:0.0" 
    #outNCFile.variables['time'].title = "Valid Time" 
    #outNCFile.variables['time'].delta_t = "000-00-00 12:00:00" 
    #outNCFile.create_variable('levids', 'c', ('level', 'chid')) 
    #outNCFile.variables['levids'].long_name = "Level Identifier" 
    #outNCFile.variables['levids'].desc = "Used to map Level with ensemble member Id" 
    outNCFile.create_variable('averageSODA', 'd', ('lat', 'lon'))
    outNCFile.variables['averageSODA'].units = "degC" 
    outNCFile.variables['lon'].assign_value(lon)
    outNCFile.variables['lat'].assign_value(lat)
    outNCFile.variables['averageSODA'].assign_value(data)
    outNCFile.close() 
  
        
def get_values():
    
    sodapath="/Volumes/HankRaid/SODA/"
   # sodapath="DATA/"
    years=[(1960+kk) for kk in range(1)]

    IDS=[(i+1) for i in range(73)]
    
    #IDS=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    No=1
    yearID=1
    log=True
    missing=["SODA_2.0.2_1958_8.cdf", "SODA_2.0.2_1958_9.cdf", "SODA_2.0.2_1959_8.cdf", "SODA_2.0.2_1959_9.cdf", "SODA_2.0.2_1981_44.cdf","SODA_2.0.2_1958_53.cdf","SODA_2.0.2_1958_63.cdf"]
    for year in years:
    
        for ID in IDS: 
            file="SODA_2.0.2_"+str(year)+"_"+str(ID)+".cdf"
            filename=sodapath+file
            
            if file not in missing:
             
                cdf_object  = IOnetcdf.read_netcdf_file(filename, log)
                temp        = cdf_object.variables["TEMP"][0,0,:,:]
                    
                 
                   
                if ID==1 and str(year)==str(years[0]):
                    lon, lat, mask_rho, temp = IOnetcdf.get_grid(cdf_object,"LON","LAT","nomask","TEMP")
                    array_size=(len(IDS),lat.shape[0],lat.shape[1])
                   
                    data=numpy.zeros((array_size),float)
                    dataNo=numpy.zeros((array_size),float)
                  
                    
                temp.resize(lat.shape)
                data[ID-1,:,:]=numpy.where(temp==-9.99e+33,data[ID-1,:,:],data[ID-1,:,:]+temp)
                dataNo[ID-1,:,:]=numpy.where(abs(temp)<50.0,dataNo[ID-1,:,:]+1.0,dataNo[ID-1,:,:]+0.0)
                
                print data[ID-1,10,10], dataNo[ID-1,10,10], temp[10,10]
                print data[ID-1,100,100], dataNo[ID-1,100,100], temp[100,100]
                
                cdf_object.close()
               
        ID=1
        yearID=yearID+1
    print 'Finished collecting all the data. Now calculate the averages:'
    
    dataFinal=numpy.zeros((lat.shape),float)
    dataFinalNo=numpy.zeros((lat.shape),float)
    
    #for k in range(len(data[1,:,1])):
    #    for kk in range(len(data[1,1,:])):
    for kkk in range(len(IDS)):
        if 55<=kkk<=73 or 0<=kkk<=18:
            
            dataFinal[:,:]=dataFinal[:,:]+(data[kkk,:,:])

            
            dataFinalNo[:,:]=dataFinalNo[:,:]+(dataNo[kkk,:,:])
    #print dataFinal[10,10], dataFinalNo[10,10], sum(dataNo[:,10,10]), sum(data[:,10,10])
    #print dataFinal[100,100], dataFinalNo[100,100], sum(dataNo[:,100,100]), sum(data[:,100,100])
    
    dataFinal=dataFinal/numpy.where(dataFinalNo!=0,dataFinalNo,1.0)
    
    print dataFinal[10,10], dataFinalNo[10,10]
    print dataFinal[100,100], dataFinalNo[100,100]
    
    plotname='average_SODA_'+str(years[0])+'s'
    return dataFinal, lat, lon, temp, plotname

def get_average():
    
    averagefile="average_SODA_1960s.cdf"
    aveNC = Nio.open_file(averagefile, "r")
    aveData = aveNC.variables["averageSODA"][:,:]
    
    return aveData

def main():
    
    """This script can be run in two ways: one to generate long-term average and write the results to a
    netCDF file using the function write_netcdf, and also to plot anomalies for a period from the
    already generated mean values, which are read from file."""
    type="Long-term average"
    #type="Anomaly"
    print "Running the script in mode: %s"%(type)
    
    """Get the data and calculate average values"""
    data, lat, lon, temp, plotname = get_values()
    
    if type=="Long-term average":
        write_netcdf(lon[1,:], lat[:,1], data, plotname)
        
    aveData = get_average()
    if type=="Anomaly":
        print data.shape, aveData.shape
        data=numpy.subtract(data,aveData)
    
    XAX=lon[1,:]
    YAX=lat[:,1]
    
    plot_soda.map_data(data,XAX,YAX,plotname)

if __name__ == "__main__":
    main()