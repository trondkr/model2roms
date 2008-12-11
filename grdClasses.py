import os, sys, datetime
from netCDF4 import Dataset
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 12, 9)
__modified__ = datetime.datetime(2008, 12, 9)
__version__  = "1.0"
__status__   = "Development"


class grdClass:

    # Note that slim and nlim are not used for mercator   
    def __init__(self, type, grdfilename,log):
        
        cdf = openNetCDF(grdfilename,log)
        
        if type=='SODA':
            print 'Generating GRD object for grid type %s'%(type)
            

    
    z_r_in = IOverticalGrid.getTypeOfVerticalStructure(vertGridSpec_IN,1,fileNameIn,log)
  
    def openNetCDF(filename, log):
        """
        Open the netCDF file and store the contents in arrays associated with variable names
        """
        try:
            cdf_file = Dataset(grdfilename,"r")
        except IOError:
            print 'Could not open file %s'%(grdfilename)
            print 'Exception caught in: openNetCDF(grdfilename, log)'
        if log is True:  
            print '\n---> Opened input file %s'%(grdfilename)
        
        return cdf_file
    
    def createObject(cdf,type):
        
        if type=='SODA':
            self.lon = cdf.variables["LON"][:]
            self.lat = cdf.variables["LAT"][:]
            self.depth = cdf.variables["DEPTH"][:]
            
                
        if type=='ROMS':
            self.lon_rho  = cdf.variables["lon_rho"][:]
            self.lat_rho  = cdf.variables["lat_rho"][:]
            self.depth    = cdf.variables["h"][:,:]
            self.mask_rho = cdf.variables["mask_rho"][:,:]       
        
        if numpy.rank(self.lon)==1:
                self.lon, self.lat = np.meshgrid(self.lon,self.lat)
            
    