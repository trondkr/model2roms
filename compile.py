
import subprocess
from datetime import datetime
import os

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 11, 11)
__modified__ = datetime(2009, 11, 11)
__version__  = "1.0"
__status__   = "Development"

def help():
    
    """
    @compile This is a simple script to call for automatic compiling of all
    fortran files necessary to run the soda2roms package. This is turned on in
    main.py with the compileAll=True
    """
    
def compileAll():
    logfile="compile.log"
    if os.path.exists(logfile): os.remove(logfile)
    log=open(logfile,'a')
    """Start the processes"""
    print "\n"
    print "Compiling interpolation.f90 to create ==> interpolation.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m interpolation interpolation.f90',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()
    log.writelines(repr(stdout_value))
    
    print "Compiling cleanArray.f90 to create ==> cl.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m cl cleanArray.f90',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))
    
    print "Compiling barotropic.f90 to create ==> barotropic.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m barotropic barotropic.f90',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))
    
    log.close()
   
    print "Compilation finished and results written to file => %s"%(logfile)
    print "\n==================================================================="

