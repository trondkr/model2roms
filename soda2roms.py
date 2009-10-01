import os, sys
from netCDF4 import Dataset
import numpy as np
import interp2D
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
import printObject
import interpolation as interp
import IOwrite
import plotData

import date
import grd
import clim2bry
import barotropic
import IOinitial
import plotData

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 8, 15)
__modified__ = datetime(2008, 8, 19)
__modified__ = datetime(2009,10, 1)
__version__  = "1.5"
__status__   = "Development"

def VerticalInterpolation(var,array1,array2,grdROMS,grdMODEL):
    
    outINDEX_ST   = (grdROMS.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
    outINDEX_U    = (grdROMS.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u,grdROMS.xi_u)
    outINDEX_V    = (grdROMS.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v,grdROMS.xi_v)

            
    if var=='salinity' or var=='temperature':
        print 'Start vertical interpolation for %s (dimensions=%s x %s)'%(var,grdROMS.xi_rho,grdROMS.eta_rho)
        outdata=np.zeros((outINDEX_ST),dtype=np.float64, order='Fortran')
    
        outdata = interp.interpolation.dovertinter(np.asarray(outdata,order='Fortran'),
                                                       np.asarray(array1,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdMODEL.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdMODEL.Nlevels),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
        
            
    if var in ['temperature','salinity']:
        #for k in range(grdROMS.Nlevels):
        #    print k
        #    plotData.contourMap(grdROMS,grdMODEL,np.squeeze(outdata[k,:,:]),k,var)
        return outdata
        
        
    if var=='vvel':
        print 'Start vertical interpolation for uvel (dimensions=%s x %s)'%(grdROMS.xi_u,grdROMS.eta_u)
        outdataU=np.zeros((outINDEX_U),dtype=np.float64,order='Fortran')
        outdataUBAR=np.zeros((outINDEX_UBAR),dtype=np.float64)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU,order='Fortran'),
                                                       np.asarray(array1,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdMODEL.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdMODEL.Nlevels),
                                                       int(grdROMS.xi_u),
                                                       int(grdROMS.eta_u),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
   
        print 'Start vertical interpolation for vvel (dimensions=%s x %s)'%(grdROMS.xi_v,grdROMS.eta_v)
        outdataV=np.zeros((outINDEX_V),dtype=np.float64,order='Fortran')
        outdataVBAR=np.zeros((outINDEX_VBAR),dtype=np.float64)
       
        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV,order='Fortran'),
                                                       np.asarray(array2,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdMODEL.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdMODEL.Nlevels),
                                                       int(grdROMS.xi_v),
                                                       int(grdROMS.eta_v),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
        
        outdataUBAR  = barotropic.velocity.ubar(np.asarray(outdataU,order='Fortran'),
                                                np.asarray(outdataUBAR,order='Fortran'),
                                                np.asarray(grdROMS.z_w,order='Fortran'),
                                                grdROMS.Nlevels,
                                                grdROMS.xi_u,
                                                grdROMS.eta_u,
                                                grdROMS.xi_rho,
                                                grdROMS.eta_rho)
        
        #plotData.contourMap(grdROMS,grdMODEL,outdataUBAR,"1","ubar")
   
        outdataVBAR  = barotropic.velocity.vbar(np.asarray(outdataV,order='Fortran'),
                                                np.asarray(outdataVBAR,order='Fortran'),
                                                np.asarray(grdROMS.z_w,order='Fortran'),
                                                grdROMS.Nlevels,
                                                grdROMS.xi_v,
                                                grdROMS.eta_v,
                                                grdROMS.xi_rho,
                                                grdROMS.eta_rho)
        
        #plotData.contourMap(grdROMS,grdMODEL,outdataVBAR,"1","vbar")
   
        return outdataU,outdataV,outdataUBAR,outdataVBAR


def HorizontalInterpolation(var,grdROMS,grdMODEL,data,show_progress):
    print 'Start %s horizontal interpolation for %s'%(grdMODEL.grdType,var)
    
    if grdMODEL.grdType=='regular':
        if var=='temperature':
            array1 = interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdMODEL,data,show_progress)
    
        if var=='salinity':
            array1 = interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdMODEL,data,show_progress)
        if var=='ssh':
            array1 = interp2D.doHorInterpolationSSHRegularGrid(var,grdROMS,grdMODEL,data)
        if var=='uvel':
            array1 = interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdMODEL,data,show_progress)
        if var=='vvel':
            array1 = interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdMODEL,data,show_progress)
    
        
    
    if grdMODEL.grdType=='irregular':
        if var=='temperature':
            interp2D.doHorInterpolationIrregularGrid(var,grdROMS,grdMODEL,data)
        if var=='salinity':
            interp2D.doHorInterpolationIrregularGrid(var,grdROMS,grdMODEL,data)
        if var=='ssh':
            interp2D.doHorInterpolationSSHIrregularGrid(var,grdROMS,grdMODEL,data)
        if var=='uvel':
            interp2D.doHorInterpolationIrregularGrid('uvel',grdROMS,grdMODEL,data)
        if var=='vvel':
            interp2D.doHorInterpolationIrregularGrid('vvel',grdROMS,grdMODEL,data)

    return array1

def rotate(grdROMS,grdMODEL,data,u,v):    
 
        """
        First rotate the values of U, V at rho points with the angle, and then interpolate
        the rho point values to U and V points and save the result
        """
        
        urot=np.zeros((int(grdMODEL.Nlevels),int(grdROMS.eta_rho),int(grdROMS.xi_rho)), np.float64)
        vrot=np.zeros((int(grdMODEL.Nlevels),int(grdROMS.eta_rho),int(grdROMS.xi_rho)), np.float64)
        
        urot, vrot = interp.interpolation.rotate(np.asarray(urot,order='Fortran'),
                                                 np.asarray(vrot,order='Fortran'),
                                                 np.asarray(u,order='Fortran'),
                                                 np.asarray(v,order='Fortran'),
                                                 np.asarray(grdROMS.angle,order='Fortran'),
                                                 int(grdROMS.xi_rho),
                                                 int(grdROMS.eta_rho),
                                                 int(grdMODEL.Nlevels))
        return urot, vrot
    
def interpolate2UV(grdROMS,grdMODEL,urot, vrot):        
        
        Zu=np.zeros((int(grdMODEL.Nlevels),int(grdROMS.eta_u),int(grdROMS.xi_u)), np.float64)
        Zv=np.zeros((int(grdMODEL.Nlevels),int(grdROMS.eta_v),int(grdROMS.xi_v)), np.float64)
        
        """
        Interpolate from RHO points to U and V points for velocities
        """

    
        Zu = interp.interpolation.rho2u(np.asarray(Zu,order='Fortran'),
                                        np.asarray(urot,order='Fortran'),
                                        int(grdROMS.xi_rho),
                                        int(grdROMS.eta_rho),
                                        int(grdMODEL.Nlevels))
        
     
        #plotData.contourMap(grdROMS,grdMODEL,Zu[0,:,:],"1",'urot')
        
        Zv = interp.interpolation.rho2v(np.asarray(Zv,order='Fortran'),
                                        np.asarray(vrot,order='Fortran'),
                                        int(grdROMS.xi_rho),
                                        int(grdROMS.eta_rho),
                                        int(grdMODEL.Nlevels))

        
        #plotData.contourMap(grdROMS,grdMODEL,Zv[0,:,:],"1",'vrot')
    
        return Zu, Zv
    
def getTime(cdf,grdROMS,grdMODEL,year,ID,type):
    
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    ref_date = date.Date()
    ref_date.day=1
    ref_date.month=1
    ref_date.year=1948
    jdref=ref_date.ToJDNumber()
    
    if type=='HYCOM':
        dateHYCOM=str(cdf.variables["Date"][0])

        hycom_date = date.Date()
        hycom_date.day=int(dateHYCOM[6:8])
        hycom_date.month=int(dateHYCOM[4:6])
        hycom_date.year=int(dateHYCOM[0:4])
        jdhycom=hycom_date.ToJDNumber()
    
        grdROMS.time=(jdhycom-jdref)
        grdROMS.reftime=jdref
       
        print '\nCurrent time of HYCOM file : %s/%s/%s'%(hycom_date.year,hycom_date.month,hycom_date.day)
        
    if type=='SODA':
        """
        Find the day and month that the SODA file respresents based on the year and ID number.
        Each SODA file represents a 5 day average, therefore we let the date we find be the first day
        of those 5 days. Thats the reason we subtract 4 below for day of month.
        """
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
    
        grdROMS.time=(jdsoda-jdref)
        grdROMS.reftime=jdref
       
        print '\nCurrent time of SODA file : %s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    if type=='SODAMONTHLY':
        """
        Find the day and month that the SODAMONTHLY file respresents based on the year and ID number.
        Each SODA file represents a 1 month average.
        """
        month=ID
        day=15
                
        soda_date = date.Date()
        soda_date.day=day
        soda_date.month=month
        soda_date.year=year
        jdsoda=soda_date.ToJDNumber()
    
        grdROMS.time=(jdsoda-jdref)
        grdROMS.reftime=jdref
       
        print '\nCurrent time of SODAMONTHLY file : %s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    
def find_subset_indices(grdMODEL,min_lat,max_lat,min_lon,max_lon):
    """
    Get the indices that covers the new grid, and enables us to only store a subset of
    the large input grid.
    """
    lat=grdMODEL.lat[:,0]
    lon=grdMODEL.lon[0,:]
 
    distances1 = []
    distances2 = []
    indices=[]
    index=0
    for point in lat:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index))
        index=index+1
        
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    
    distances1 = []
    distances2 = []
    index=0
   
    for point in lon:
        s1 = max_lon-point # (vector subtract)
        s2 = min_lon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index))
        index=index+1
       
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    
    """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
      
    grdMODEL.minJ=indices[1][2]
    grdMODEL.maxJ=indices[0][2]
    grdMODEL.minI=indices[3][2]
    grdMODEL.maxI=indices[2][2]
    
  
def convertMODEL2ROMS(years,IDS,climName,initName,dataPath,romsgridpath,vars,show_progress,type):

    if type=='SODA':
        fileNameIn=dataPath+'SODA_2.0.2_'+str(years[0])+'_'+str(IDS[0])+'.cdf'
    if type=='SODAMONTHLY':
        fileNameIn=dataPath+'SODA_2.0.2_'+str(years[0])+'0'+str(IDS[0])+'.cdf'
        
    if type=='HYCOM':
       fileNameIn=dataPath+'archv.2003_307_00_3zt.nc'
       
    """
    First time in loop, get the essential old grid information
    MODEL data already at Z-levels. No need to interpolate to fixed depths,
    but we use the one we have
    """
    
    grdMODEL = grd.grdClass(fileNameIn,type)
    grdROMS = grd.grdClass(romsgridpath,"ROMS")
    
    """Now we want to subset the data to avoid storing more information than we need.
    We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes"""
    if type=='SODA' or type=='SODAMONTHLY':
        find_subset_indices(grdMODEL,min_lat=30, max_lat=90, min_lon=0, max_lon=360)
    if type=='HYCOM':
        #GOM
        grdMODEL.minJ=1900
        grdMODEL.maxJ=2400
        grdMODEL.minI=2575
        grdMODEL.maxI=2875
     

    grdMODEL.lat=grdMODEL.lat[grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI]
    grdMODEL.lon=grdMODEL.lon[grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI]
  
    print "\n---> Selected area in output file spans from (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)"%(grdROMS.lon_rho.min(),grdROMS.lat_rho.min(),grdROMS.lon_rho.max(),grdROMS.lat_rho.max())
    print "---> Selected area in input file spans from  (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)\n"%(grdMODEL.lon.min(),grdMODEL.lat.min(),grdMODEL.lon.max(),grdMODEL.lat.max())
    
    print '\n---> Finished initializing'
    print '\n--------------------------\n'
    
    time=0
   
    
    for year in years:
        
        firstRun = True ;
        
        for ID in IDS:
           
            if type=='SODA':
                file="SODA_2.0.2_"+str(year)+"_"+str(ID)+".cdf"
                filename=dataPath+file
                variableNames=['TEMP','SALT','SSH','U','V']
                
            if type=='SODAMONTHLY':
                if ID <  10: filename=dataPath+'SODA_2.0.2_'+str(years[0])+'0'+str(ID)+'.cdf'
                if ID >= 10: filename=dataPath+'SODA_2.0.2_'+str(years[0])+str(ID)+'.cdf'
                variableNames=['temp','salt','ssh','u','v']
              
            if type=='HYCOM':
                filename=dataPath+'archv.2003_307_00_3zt.nc'
                variableNames=['temperature','salinity','SSH','U','V'] # NATHAN FIXME; give correct name of hycom variables
                
            """Now open the input file"""    
            cdf = Dataset(filename)
            
            getTime(cdf,grdROMS,grdMODEL,year,ID,type)
            
            """Each MODEL file consist only of one time step. Get the subset data selected, and
            store that time step in a new array:"""
            
            if firstRun is True:
                firstRun = False
   
                indexROMS_S_ST = (grdROMS.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
                indexROMS_SSH  = (grdROMS.eta_rho,grdROMS.xi_rho)
                
                indexROMS_S_U = (grdROMS.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
                indexROMS_S_V = (grdROMS.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
                indexROMS_UBAR = (grdROMS.eta_u,grdROMS.xi_u)
                indexROMS_VBAR = (grdROMS.eta_v,grdROMS.xi_v)
        
            """
            All variables for all time are now stored in arrays. Now, start the interpolation to the
            new grid for all variables and then finally write results to file.
            """
            
            for var in vars:
                
                if var=='temperature':
                    STdata=np.zeros((indexROMS_S_ST),dtype=np.float64)
                    if type=="SODA": data = np.array(cdf.variables[str(variableNames[0])][0,:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    if type=="SODAMONTHLY": data = np.array(cdf.variables[str(variableNames[0])][:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    
                    if time==0:
                        tmp=np.squeeze(data[0,:,:])
                        grdMODEL.mask = np.zeros((grdMODEL.lon.shape),dtype=np.float64)
                        grdMODEL.mask[:,:] = np.where(tmp==grdROMS.fill_value,1,0)
                        
                if var=='salinity':
                    STdata = np.zeros((indexROMS_S_ST),dtype=np.float64)
                    if type=="SODA": data  = np.array(cdf.variables[str(variableNames[1])][0,:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    if type=="SODAMONTHLY": data  = np.array(cdf.variables[str(variableNames[1])][:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                
                if var=='ssh':
                    SSHdata = np.zeros((indexROMS_SSH),dtype=np.float64)
                    if type=="SODA": data  = np.array(cdf.variables[str(variableNames[2])][0,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    if type=="SODAMONTHLY": data  = np.array(cdf.variables[str(variableNames[2])][grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                
                if var=='uvel':
                    Udata = np.zeros((indexROMS_S_U),dtype=np.float64)
                    if type=="SODA": data  = np.array(cdf.variables[str(variableNames[3])][0,:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    if type=="SODAMONTHLY": data  = np.array(cdf.variables[str(variableNames[3])][:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    
                if var=='vvel':
                    Vdata = np.zeros((indexROMS_S_V),dtype=np.float64)
                    if type=="SODA": data  = np.array(cdf.variables[str(variableNames[4])][0,:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    if type=="SODAMONTHLY": data  = np.array(cdf.variables[str(variableNames[4])][:,grdMODEL.minJ:grdMODEL.maxJ,grdMODEL.minI:grdMODEL.maxI])
                    
                array1 = HorizontalInterpolation(var,grdROMS,grdMODEL,data,show_progress)
                
                if var in ['temperature','salinity']:    
                    STdata = VerticalInterpolation(var,array1,array1,grdROMS,grdMODEL)
                
                    IOwrite.writeClimFile(grdROMS,time,climName,var,STdata)
                    if time==0 and grdROMS.write_init is True:
                        IOinitial.createInitFile(grdROMS,time,initName,var,STdata)
                        
                if var=='ssh':
                    SSHdata=array1[0,:,:]
                    IOwrite.writeClimFile(grdROMS,time,climName,var,SSHdata)
                    if time==0:
                        IOinitial.createInitFile(grdROMS,time,initName,var,SSHdata)
                        
                if var=='uvel':
                    array2=array1
                    
                if var=='vvel':
        
                    UBARdata = np.zeros((indexROMS_UBAR),dtype=np.float64)
                    VBARdata = np.zeros((indexROMS_VBAR),dtype=np.float64)
                
                    urot,vrot = rotate(grdROMS,grdMODEL,data,array2,array1)
                    
                    u,v = interpolate2UV(grdROMS,grdMODEL,urot,vrot)
                   
                    Udata,Vdata,UBARdata,VBARdata = VerticalInterpolation(var,u,v,grdROMS,grdMODEL)
                
                if var=='vvel':
                    IOwrite.writeClimFile(grdROMS,time,climName,var,Udata,Vdata,UBARdata,VBARdata)
                    if time==0:
                        IOinitial.createInitFile(grdROMS,time,initName,var,Udata,Vdata,UBARdata,VBARdata)
                  
                
            cdf.close()    
            if show_progress is True:
                from progressBar import progressBar
                # find unicode characters here: http://en.wikipedia.org/wiki/List_of_Unicode_characters#Block_elements
                empty  =u'\u25FD'
                filled =u'\u25FE'
        
                progress = progressBar(color='green',width=24, block=filled.encode('UTF-8'), empty=empty.encode('UTF-8'))
                message='Finished conversions for time %s'%(grdROMS.message)
                progress.render(100,message)
            else:
                message='Finished conversions for time %s'%(grdROMS.message)
                print message
            time+=1
           
          
            
    
