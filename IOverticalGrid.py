
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

    
def get_dims(var):
    dims=var.shape
    return dims

def calculate_z_w(self):
    
    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[eta_rho, xi_rho, s]
     
    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008
    """
    
    h =self.depth
    
    hc=np.min(h,self.Tcline)
    
    if self.theta_s != 0.0:
        cff1=1.0/np.sinh(self.theta_s)
        cff2=0.5/np.tanh(0.5*self.theta_s)

    cff=1.0/self.Nlevels
    
    sc_w=np.zeros((self.Nlevels+1),float)
    Cs_w=np.zeros((self.Nlevels+1),float)
    
    sc_w[0]=-1.0
    Cs_w[0]=-1.0
    
    for k in xrange(self.Nlevels):
        sc_w[k+1]=cff*(float(k-self.Nlevels))
        if self.theta_s != 0.0:
            Cs_w[k+1]=(1.0-self.theta_b)*cff1*np.sinh(self.theta_s*sc_w[k+1])+self.theta_b*(cff2*np.tanh(self.theta_s*(sc_w[k+1]+0.5))- 0.5)
        else:
            Cs_w[k+1]=sc_w[k+1]
    
    """
    TODO: FIXME hardcode variables should be read from file
    """
    zeta=0.0      
    dim=get_dims(h)
    eta_rho=dim[0]
    xi_rho=dim[1]
    z_w=np.zeros((len(sc_w),eta_rho,xi_rho),float)
   
    for k in xrange(len(sc_w)):
        z_w[k,:,:] = np.multiply(zeta,(1+sc_w[k]))+hc*(sc_w[k])+(np.subtract(h,hc))*Cs_w[k]
      
    self.z_w = z_w
    self.Cs_w=Cs_w
    self.s_w=sc_w
    

def calculate_z_r(self):
    
    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[eta_rho, xi_rho, s]
     
    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008
    """
    
    h =self.depth
    
    hc=np.min(h,self.Tcline)
    
    if self.theta_s != 0.0:
        cff1=1.0/np.sinh(self.theta_s)
        cff2=0.5/np.tanh(0.5*self.theta_s)

    cff=1.0/self.Nlevels
    
    sc_r=np.zeros((self.Nlevels),float)
    
            
    Cs_r=np.zeros((self.Nlevels),float)
    
    for k in xrange(self.Nlevels):
        sc_r[k]=cff*(float(k-self.Nlevels)-0.5)
        if self.theta_s != 0.0:
            Cs_r[k]=(1.0-self.theta_b)*cff1*np.sinh(self.theta_s*sc_r[k])+self.theta_b*(cff2*np.tanh(self.theta_s*(sc_r[k]+0.5))- 0.5)
        else:
            Cs_r[k]=sc_r[k]
    
    """
    TODO: FIXME hardcode variables should be read from file
    """
    zeta=0.0      
    dim=get_dims(h)
    eta_rho=dim[0]
    xi_rho=dim[1]
    z_r=np.zeros((len(sc_r),eta_rho,xi_rho),float)
   
    for k in range(len(sc_r)):
        z_r[k,:,:] = np.multiply(zeta,(1+sc_r[k]))+hc*(sc_r[k])+(np.subtract(h,hc))*Cs_r[k]
      
    self.z_r = np.flipud(z_r)
    self.Cs_rho=Cs_r
    self.hc=hc
    self.s_rho=sc_r
