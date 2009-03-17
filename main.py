
import os, sys, time
from datetime import datetime
import soda2roms, IOstation
import clim2bry
import grd

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009, 2,13)
__version__  = "1.1"
__status__   = "Development"

def help():
    """
    This program is run by typing: python main;py in the command window.
    """
    
def main():
    print '\n--------------------------\n'
    print '\n---> Initializing\n'
    print 'Started ' + time.ctime(time.time())
    
    climName='test.nc'
    initName='test_init.nc'
    bryName='test_bry.nc'
    
    sodapath="/Volumes/HankRaid/SODA/"
    romsgridpath="/Users/trond/Projects/arcwarm/nordic/AA_10km_grid.nc"
    romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd.nc"
    start_year=1960
    end_year=1961
    start_day_in_start_year=10
    end_day_in_end_year=65
    
    years=[(int(start_year)+kk) for kk in range(int(end_year)-int(start_year))]
    
    IDS=[(0+i+1) for i in range(73)]
    
    soda2roms.convertSODA2ROMS(years,IDS,climName,initName,sodapath,romsgridpath)
    
    
    grdROMS = grd.grdClass(romsgridpath,"ROMS")
    
    clim2bry.writeBry(grdROMS,'1960',bryName,climName)   

    # GB, NovaScotia, Grand Bank, Nuuk, Iceland, NS, Lofoten, BS
   
        
    lonlist=[-66.5]#,-66.40,-50.43,-54.38, 21.51,  1.53, 13.38, 32.75]
    latlist=[ 41.5]#, 43.41, 44.50, 64.72, 63.35, 58.36, 67.50, 71.77]
    #IOstation.getStationData(years,IDS,outfilename,sodapath,latlist,lonlist)

    print 'Finished ' + time.ctime(time.time())
    

if __name__ == "__main__":
    #import profile
    #
    #try:
    #    import psyco
    #    psyco.log()
    #    psyco.profile(0.2)
    #
    #except ImportError:
    #    pass
    #profile.run('main()')
    
    main()
