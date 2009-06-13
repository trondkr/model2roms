
import os, sys, time
from datetime import datetime
import soda2roms, IOstation
import clim2bry, DecimateGrid
import grd
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009, 4,22)
__version__  = "1.1"
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
    convert2Grid=True
    decimateGrid=False 
    
    sodapath="/Volumes/HankRaid/SODA/"
    hycompath="/Users/trond/Projects/arcwarm/SODA/HYCOM/"
    
    #romsgridpath="/Users/trond/Projects/arcwarm/nordic/AA_10km_grid.nc"
    romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd.nc"
    #romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd_030208.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/nordic/imr_nordic_4km.nc"
    #romsgridpath="/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc"
    #romsgridpath="/Users/trond/Projects/Nathan/NoMed47_GRID_Global.nc"
    #romsgridpath='/Users/trond/Projects/Nathan/GOM_GRID_Global.nc'
    start_year      =1990
    end_year        =2000
    start_julianday =0
    end_julianday   =365
    
    """Set the input data MODEL type: Current options are SODA or HYCOM"""
    type='HYCOM' 
    type='SODA'
    
    vars=['temperature','salinity','ssh','uvel','vvel']
    vars=['temperature']
    
    
    start_day_in_start_year=np.round(start_julianday/5.0)
    end_day_in_end_year=round(end_julianday/5.0)
    
    years=[(int(start_year)+kk) for kk in range(int(end_year)-int(start_year))]
    
    loop=int(end_day_in_end_year)-int(start_day_in_start_year)
   
    if int(start_day_in_start_year)==int(end_day_in_end_year):
        IDS=[int(start_day_in_start_year)]
    else:
        IDS=[(i+1+int(start_day_in_start_year)) for i in range(loop)]
    
    grdROMS = grd.grdClass(romsgridpath,"ROMS")
    
    if convert2Grid==True:
        
        showInfo(vars,romsgridpath,climName,initName,bryName,start_year,end_year,start_julianday,end_julianday)
        
        if type=='SODA':
            soda2roms.convertMODEL2ROMS(years,IDS,climName,initName,sodapath,romsgridpath,vars,show_progress,type)
        elif type=='HYCOM':
            soda2roms.convertMODEL2ROMS([1],[1,2,3],climName,initName,hycompath,romsgridpath,vars,show_progress,type)
    
        clim2bry.writeBry(grdROMS,'1960',bryName,climName)
    
    
    
    if decimateGrid==True:
        DecimateGrid.createGrid(grdROMS,'/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc',2)
     

    if write_stations is True:
        print "Running in station mode"
        # GB, NovaScotia, Grand Bank, Nuuk, Iceland, NS, Lofoten, BS     
        lonlist=[-66.5233]#,-66.40,-50.43,-54.38, 21.51,  1.53, 13.38, 32.75]
        latlist=[ 41.5423]#, 43.41, 44.50, 64.72, 63.35, 58.36, 67.50, 71.77]
        
        IOstation.getStationData(years,IDS,sodapath,latlist,lonlist)

    print 'Finished ' + time.ctime(time.time())
    

if __name__ == "__main__":
    import psyco
    try:
        import psyco
        psyco.log()
        psyco.full(memory=100)
        psyco.profile(0.05, memory=100)
        psyco.profile(0.2)
    except ImportError:
        pass

    main()
