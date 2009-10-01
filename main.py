
import os, sys, time
from datetime import datetime
import soda2roms, IOstation
import clim2bry, DecimateGrid
import grd
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009,10, 1)
__version__  = "1.2"
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
    
    climName='gom_clim_1990.nc'
    initName='gom_init_1990.nc'
    bryName='gom_bry_1990.nc'
    
    """Set show_progress to "False" if you do not want to see the progress
    indicator for horizontal interpolation. This requires the two modules:
    terminal.py and progressBar.py"""
    show_progress=True
    write_stations=False
    createForcing=True
    decimateGrid=False 
    
    sodapath="/Volumes/HankRaid/SODA/"
    sodapath="/Volumes/MacintoshHD2/Datasets/SODA/"
    sodapath="/Volumes/MacintoshHD2/Datasets/SODAMonthly/"
    hycompath="/Users/trond/Projects/arcwarm/SODA/HYCOM/"
    
    #romsgridpath="/Users/trond/Projects/arcwarm/nordic/AA_10km_grid.nc"
    romsgridpath="/Volumes/HankRaid/ROMS/GoM/grid/gom_grd.nc"
    #romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd_030208.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/nordic/imr_nordic_4km.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc"
    #romsgridpath="/Users/trond/Projects/Nathan/NoMed47_GRID_Global.nc"
    #romsgridpath='/Users/trond/Projects/Nathan/GOM_GRID_Global.nc'
    start_year      =1959
    end_year        =1964
    start_julianday =0
    end_julianday   =365
    
    """Set the input data MODEL type: Current options are SODA or HYCOM"""
    type='HYCOM' 
    type='SODA'
    type='SODAMONTHLY'
    
    vars=['temperature','salinity','ssh','uvel','vvel']
    #vars=['temperature']
    
    if type=='SODA': aveDays=5.0
    if type=='SODAMONTHLY': aveDays=30.0
    
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
    
        clim2bry.writeBry(grdROMS,'1960',bryName,climName)
    
    
    
    if decimateGrid==True:
        DecimateGrid.createGrid(grdROMS,'/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc',2)
     

    if write_stations is True:
        print "Running in station mode"
        stationNames=['NorthSea','Iceland','EastandWestGreenland','Lofoten', 'Georges Bank']
        lonlist=[ 2.4301, -22.6001, -47.0801,  13.3801, -67.2001]
        latlist=[54.5601, 63.7010,  60.4201,  67.5001,  41.6423]
    
    
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
