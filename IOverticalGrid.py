
from datetime import datetime
import numpy as np


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 8, 15)
__modified__ = datetime(2008, 8, 19)
__modified__ = datetime(2008, 12, 9)
__version__  = "1.1"
__status__   = "Development"


def get_z_levels(self):
    """
    Get a list of all the variables contained in netCDF file "filename"
    """
    if self.type=='SODA':
            self.z_r=-self.depth
    if len(self.z_r)==0:
        print "No depth matrix found in file %s"%(self.selffilename)


def calculate_z_w(self):
    
    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[eta_rho, xi_rho, s]
     
    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008
    """
    sc_w=np.zeros((self.Nlevels+1),np.float64)
    Cs_w=np.zeros((self.Nlevels+1),np.float64)
    z_w=np.zeros((len(sc_w),self.eta_rho,self.xi_rho),np.float64)
    
    h =self.depth
    hc=np.min(h,self.Tcline)
    cff=1.0/float(self.Nlevels+1)
    
    for k in xrange(self.Nlevels+1):
        sc_w[k]=-1.0+(self.Nlevels-k)*cff
    
        Ptheta=np.sinh(self.theta_s*sc_w[k])/np.sinh(self.theta_s)
        Rtheta=np.tanh(self.theta_s*(sc_w[k]+0.5))/(2.0*np.tanh(0.5*self.theta_s))-0.5
    
        if self.theta_s != 0.0:
            Cs_w[k]=(1-self.theta_b)*Ptheta+self.theta_b*Rtheta
        else:
            Cs_w[k]=sc_w[k]
            
   
    """
    TODO: FIXME hardcode variables should be read from file
    """
    zeta=None
    
    for k in xrange(len(sc_w)):
        scmCshc=(np.subtract(sc_w[k],Cs_w[k]))*hc
        z_w[k,:,:] = scmCshc + np.multiply(Cs_w[k],h)
        
        if zeta != None:
            dd = np.divide(zeta,h)
            z_w[k,:,:] = z_w[k,:,:] + scmCshc*dd
            
    self.z_w = z_w
    self.Cs_w=Cs_w
    self.hc=hc
    self.s_w=sc_w
    
    #for l in xrange(10,400,20):
    #    for ll in xrange(10,300,20):
    #        print z_w[:,ll,l],self.depth[ll,l]
  
def calculate_z_r(self):
    
    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[s,eta_rho,xi_rho]
     
    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008, 11.03.2009
    """
    sc_r=np.zeros((self.Nlevels),np.float64)
    Cs_r=np.zeros((self.Nlevels),np.float64)
    z_r=np.zeros((len(sc_r),self.eta_rho,self.xi_rho),np.float64)
    
    h =self.depth
    hc=np.min(h,self.Tcline)
    cff=1.0/float(self.Nlevels)
            
    for k in xrange(self.Nlevels):
        sc_r[k]=-1.0+(self.Nlevels-k-0.5)*cff
    
        Ptheta=np.sinh(self.theta_s*sc_r[k])/np.sinh(self.theta_s)
        Rtheta=np.tanh(self.theta_s*(sc_r[k]+0.5))/(2.0*np.tanh(0.5*self.theta_s))-0.5
    
        if self.theta_s != 0.0:
            Cs_r[k]=(1-self.theta_b)*Ptheta+self.theta_b*Rtheta
        else:
            Cs_r[k]=sc_r[k]

    """
    TODO: FIXME hardcode variables should be read from file
    """
    zeta=None
    
    for k in xrange(len(sc_r)):
        scmCshc=(np.subtract(sc_r[k],Cs_r[k]))*hc
        z_r[k,:,:] = scmCshc + np.multiply(Cs_r[k],h)
        
        if zeta != None:
            dd = np.divide(zeta,h)
            z_r[k,:,:] = z_r[k,:,:] + scmCshc*dd
  
    self.z_r = z_r
    self.Cs_rho=Cs_r
    self.hc=hc
    self.s_rho=sc_r
    