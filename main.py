
import os, sys, time
from datetime import datetime
import soda2roms, IOstation
import clim2bry, DecimateGrid
import grd
import numpy as np
import ncarRiver
#import ncarSSS


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009,11, 9)
__version__  = "1.3"
__status__   = "Development"

def help():
    """
    This program is run by typing: python main.py in the command window.
    """
    
def showInfo(vars,romsgridpath,climName,initName,bryName,start_year,end_year,start_julianday,end_julianday):
    print 'Conversions run from day %s in %s to day %s in year %s'%(start_julianday,start_year,end_julianday,end_year)
    print 'The following variables will be converted:'
    for var in vars:
        print '---> %s'%(var)
    print '\nOutput grid file is: %s'%(romsgridpath)
    print '\n--------------------------\n'
    print '\n---> Initializing'
    
    
def main():
    print '\n--------------------------\n'
    print 'Started ' + time.ctime(time.time())
    
    """ EDIT ==================================================================="""
    
    """Set show_progress to "False" if you do not want to see the progress
    indicator for horizontal interpolation. This requires the two modules:
    terminal.py and progressBar.py"""
    show_progress=True
    """Set compileAll to True if you want automatic re-compilation of all the
    fortran files necessary to run soda2roms. You need to edit compile.py for this"""
    compileAll=False
    """Extract time-series of data for given longitude/latitude""" 
    extractStations=False
    """Create the bry, init, and clim files for a given grid and input data"""
    createForcing=True
    """Create a smaller resolution grid based on your original. Decimates every second for
    each time run"""
    decimateGrid=False
    """Create river runoff file based on NCEP-NCAR CORE data"""
    createRiverRunoff=False
    
    """Define the paths to the CORE data (only river file), and the SODA/HYCOM data"""
    corepath="/Users/trond/Projects/arcwarm/CORE/"
    sodapath="/Volumes/HankRaid/SODA/"
    sodapath="/Volumes/MacintoshHD2/Datasets/SODA/"
    sodapath="/Volumes/MacintoshHD2/Datasets/SODAMonthly/"
    sodapath="/Users/trond/Projects/arcwarm/SODAmonthly/"
    hycompath="/Users/trond/Projects/arcwarm/SODA/HYCOM/"
    
    romsgridpath="/Volumes/HankRaid/ROMS/GoM/grid/gom_grd.nc"
    romsgridpath="/Users/trond/Projects/Roms/GOMfull/Inputs/gom_grd.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/nordic/imr_nordic_4km.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc"
    
    start_year      =1960
    end_year        =1961
    start_julianday =0
    end_julianday   =365
    
    """Name of output files for CLIM, BRY, and INIT files"""
    climName='gom_clim_'+str(start_year)+'.nc'
    initName='gom_init_'+str(start_year)+'.nc'
    bryName='gom_bry_'+str(start_year)+'.nc'
    
    """Set the input data MODEL type: Current options are SODA or HYCOM"""
    type='HYCOM' 
    type='SODA'
    type='SODAMONTHLY'
    
    """Define what variables to include in the forcing files"""
    vars=['temperature','salinity','ssh','uvel','vvel']
    #vars=['uvel','vvel']
    
    """5 day or 30 day average files for SODA"""
    if type=='SODA': aveDays=5.0
    if type=='SODAMONTHLY': aveDays=30.0
    
    """Define a set of longitude/latitude positions with names to extract into
    station files (using extractStations)"""
    stationNames=['NorthSea','Iceland','EastandWestGreenland','Lofoten', 'Georges Bank']
    lonlist=[ 2.4301, -22.6001, -47.0801,  13.3801, -67.2001]
    latlist=[54.5601, 63.7010,  60.4201,  67.5001,  41.6423]
    
        
    """" NO EDIT BELOW ========================================================="""
    if compileAll is True:
        import compile
        compile.compileAll()
        
    start_day_in_start_year=np.round(start_julianday/aveDays)
    end_day_in_end_year=round(end_julianday/aveDays)
    
    years=[(int(start_year)+kk) for kk in range(int(end_year)-int(start_year))]
    
    loop=int(end_day_in_end_year)-int(start_day_in_start_year)
   
    if int(start_day_in_start_year)==int(end_day_in_end_year):
        IDS=[int(start_day_in_start_year)]
    else:
        IDS=[(i+1+int(start_day_in_start_year)) for i in range(loop)]
    
    grdROMS = grd.grdClass(romsgridpath,"ROMS")

    if createForcing==True:
        
        showInfo(vars,romsgridpath,climName,initName,bryName,start_year,end_year,start_julianday,end_julianday)
        
        if type=='SODA' or type=='SODAMONTHLY':
            soda2roms.convertMODEL2ROMS(years,IDS,climName,initName,sodapath,romsgridpath,vars,show_progress,type)
        elif type=='HYCOM':
            soda2roms.convertMODEL2ROMS([1],[1,2,3],climName,initName,hycompath,romsgridpath,vars,show_progress,type)
    
        clim2bry.writeBry(grdROMS,start_year,bryName,climName)
    
    if createRiverRunoff is True:
        ncarRiver.runoff(grdROMS,corepath)
    
    
    if decimateGrid==True:
        DecimateGrid.createGrid(grdROMS,'/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc',2)
     

    if extractStations is True:
        print "Running in station mode and extracting pre-defined station locations"
        IOstation.getStationData(years,IDS,sodapath,latlist,lonlist,stationNames)

    print 'Finished ' + time.ctime(time.time())
    

if __name__ == "__main__":
    
    try:
        import psyco
        psyco.log()
        psyco.full(memory=100)
        psyco.profile(0.05, memory=100)
        psyco.profile(0.2)
    except ImportError:
        pass

    main()
