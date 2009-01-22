
import os, sys, datetime
import numpy

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2008, 8, 19)
__modified__ = datetime.datetime(2009, 1, 22)
__version__  = "1.0"
__status__   = "Development"

def get_station(grdSODA,st_lon,st_lat):
    """
    This is a function that takes longitude and latitude as
    decimal input, and returns the index values closest to
    the longitude and latitude. This is an iterative process for finding the best
    index pair. Trond Kristiansen, 11.03.2008
    """
    longitude, latitude, mask_rho, depth = IOnetcdf.get_grid(cdf_object,"LON","LAT","nomask","DEPTH")
    #Find closest latitude column 
    tmp2=numpy.abs(subtract(latitude,st_lat))
    #Sort the column to find the index that closest matches station
    delta_lat = tmp2.argsort(axis=1)

    # Find closest longitude index in latitude column
    tmp1=numpy.abs(subtract(longitude[:,delta_lat[0,0]],st_lon))
    # Sort the result to find the best longitude in that column
    delta_lon = tmp1.argsort()
    tmp2=numpy.abs(subtract(latitude[delta_lon[0],:],st_lat))
    # Repeat process to find the best latitude in the best longitude column
    delta_lat = tmp2.argsort(axis=0)
    tmp1=numpy.abs(subtract(longitude[:,delta_lat[0]],st_lon))

    delta_lon = tmp1.argsort()
    
    if (stdlog is True):
        print '' 
        print '=====get_station======'
        print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon,st_lat)
        print 'Result ===>'
        print 'Found index pair [%s,%s] in gridfile'%(delta_lon[0],delta_lat[0])
        print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[delta_lon[0],delta_lat[0]],latitude[delta_lon[0],delta_lat[0]])
        print '======================'
        print ''
        
    return delta_lon[0], delta_lat[0]

