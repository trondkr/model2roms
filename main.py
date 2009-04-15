
import os, sys, time
from datetime import datetime
import soda2roms, IOstation
import clim2bry
import grd
import numpy as np

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009, 2,13)
__version__  = "1.1"
__status__   = "Development"

def help():
    """
    This program is run by typing: python main.py in the command window.
    """
    
def main():
    print '\n--------------------------\n'
    print 'Started ' + time.ctime(time.time())
    
    climName='test.nc'
    initName='test_init.nc'
    bryName='test_bry.nc'
    
    write_clim=True
    write_bry=False
    write_init=True
    write_stations=False
    
    """Set show_progress to "False" if you do not want to see the progress
    indicator for horizontal interpolation. This requires the two modules:
    terminal.py and progressBar.py"""
    show_progress=True
    
    sodapath="/Volumes/HankRaid/SODA/"
    romsgridpath="/Users/trond/Projects/arcwarm/nordic/AA_10km_grid.nc"
    romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd.nc"
  #  romsgridpath="/Users/trond/Projects/arcwarm/nordic/imr_nordic_4km.nc"
    start_year      =1960
    end_year        =1961
    start_julianday =1
    end_julianday   =50
    
    vars=['temperature','salinity','ssh','uvel','vvel']
    vars=['temperature']
    print 'Conversions run from day %s in %s to day %s in year %s'%(start_julianday,start_year,end_julianday,end_year)
    print 'The following variables will be converted:'
    for var in vars:
        print '---> %s'%(var)
    print '\nOutput grid file is: %s'%(romsgridpath)
    print '\nThe following functions are run:'
    if write_clim is True: print '---> write clim (%s)'%(climName)
    if write_init is True: print '---> write init (%s)'%(initName)
    if write_bry is True: print '---> write bry (%s)'%(bryName)
    if write_stations is True: print '---> write stations'
    print '\n--------------------------\n'
    print '\n---> Initializing'
    
    start_day_in_start_year=np.round(start_julianday/5.0)
    end_day_in_end_year=round(end_julianday/5.0)
    
    years=[(int(start_year)+kk) for kk in range(int(end_year)-int(start_year))]

    IDS=[(i+1+int(start_day_in_start_year)) for i in range(int(end_day_in_end_year))]
    
    soda2roms.convertSODA2ROMS(years,IDS,climName,initName,sodapath,romsgridpath,vars,show_progress)
    
    if write_bry is True:
        grdROMS = grd.grdClass(romsgridpath,"ROMS")
        clim2bry.writeBry(grdROMS,'1960',bryName,climName)   

    if write_stations is True:
        # GB, NovaScotia, Grand Bank, Nuuk, Iceland, NS, Lofoten, BS     
        lonlist=[-66.5]#,-66.40,-50.43,-54.38, 21.51,  1.53, 13.38, 32.75]
        latlist=[ 41.5]#, 43.41, 44.50, 64.72, 63.35, 58.36, 67.50, 71.77]
        IOstation.getStationData(years,IDS,outfilename,sodapath,latlist,lonlist)

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
