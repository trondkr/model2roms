
import os, sys, time
from datetime import datetime
import soda2roms, IOstation


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 1,30)
__modified__ = datetime(2009, 1,30)
__version__  = "1.0"
__status__   = "Development"


def main():
    print '\n--------------------------\n'
    print '\n---> Initializing\n'
    print 'Started ' + time.ctime(time.time())
    
    outfilename='test.nc'
    sodapath="/Volumes/HankRaid/SODA/"
    romsgridpath="/Users/trond/ROMS/GoM/grid/gom_grd.nc"
    start_year=1965
    end_year=1976
    start_day_in_start_year=10
    end_day_in_end_year=65
    
    years=[(int(start_year)+kk) for kk in range(int(end_year)-int(start_year))]
    
    IDS=[(0+i+1) for i in range(73)]
    
    #soda2roms.convertSODA2ROMS(years,IDS,outfilename,sodapath,romsgridpath)
    
    lonlist=[300]
    latlist=[42]
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
