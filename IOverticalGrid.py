
from datetime import datetime
import numpy as np


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 8, 15)
__modified__ = datetime(2009, 11, 11)
__version__  = "1.3"
__status__   = "Development"


def get_z_levels(self):
    """
    Get a list of all the variables contained in netCDF file "filename"
    """
    self.z_r=-self.depth

    if len(self.z_r)==0:
        print "No depth matrix found in file %s"%(self.selffilename)


def calculate_z_w(self):

    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[eta_rho, xi_rho, s]

    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008, 11.11.2009
    """
    sc_w=np.zeros((self.Nlevels+1),np.float64)
    Cs_w=np.zeros((self.Nlevels+1),np.float64)
    z_w=np.zeros((len(sc_w),self.eta_rho,self.xi_rho),np.float64)

    h =self.depth
    hc=h.min() #np.min(np.min(h),self.Tcline)
    cff=1.0/float(self.Nlevels+1)

    for k in xrange(self.Nlevels+1):
        sc_w[k]=-1.0+(self.Nlevels-k)*cff

        if self.vstretching==1:
            """
            Original vertical strectching function, Song and Haidvogel (1994).
            see: https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
            C(s) = (1 - b) * [SINH(s * a) / SINH(a)] +

             b * [-0.5 + 0.5 * TANH(a * (s + 0.5)) / TANH(0.5 * a)]
            """
            Ptheta=np.sinh(self.theta_s*sc_w[k])/np.sinh(self.theta_s)
            Rtheta=np.tanh(self.theta_s*(sc_w[k]+0.5))/(2.0*np.tanh(0.5*self.theta_s))-0.5

            if self.theta_s != 0.0:
                Cs_w[k]=(1-self.theta_b)*Ptheta+self.theta_b*Rtheta
            else:
                Cs_w[k]=sc_w[k]
            
            if k==0:
                Cs_w[k]=0
            
        if self.vstretching==2:
            """
            A. Shchepetkin new vertical stretching function.
            C(s) = [1.0 - COSH(theta_s * s)] / [COSH(theta_s) - 1.0]
            Adapted from set_scoord.F (ROMS source code)
            """
            Aweight=1.0
            Bweight=1.0

            if self.theta_s > 0.0:
                Csur=(1.0 - np.cosh(self.theta_s * sc_w[k]))/(np.cosh(self.theta_s) - 1.0)

                if self.theta_b > 0.0:

                    Cbot=np.sinh(self.theta_b)*(sc_w[k] + 1.0)/np.sinh(self.theta_b) - 1.0

                    Cweight=(sc_w[k]+1.0)**Aweight*(1.0 + (Aweight/Bweight)*(1.0-(sc_w[k]+1.0)**Bweight))

                    Cs_w[k]=Cweight*Csur +(1.0-Cweight)*Cbot
                else:
                    Cs_w[k]=Csur
            else:
                Cs_w[k]=sc_w


    zeta=None

    for k in xrange(len(sc_w)):
        if self.vstretching==1:
            z_w[k,:,:] =np.multiply(sc_w[k],hc) + np.subtract(h,hc)*Cs_w[k]
           
            
        if self.vstretching==2:
            z_w[k,:,:] =np.multiply(sc_w[k],hc) + np.subtract(h,hc)*Cs_w[k]
    
    self.z_w = z_w * self.mask_rho
    self.Cs_w=Cs_w
    self.s_w=sc_w

def calculate_z_r(self):

    """
    Function that estimates the matrix that converts sigma
    depth values to depth in meters. This matrix is time dependent
    (zeta varies), and also position dependent since the bottom matrix varies.
    Results are stored in array z[s,eta_rho,xi_rho]

    Trond Kristiansen, 20.01.2008, 03.12.2008, 09.12.2008, 11.03.2009, 11.11.2009
    """
    sc_r=np.zeros((self.Nlevels),np.float64)
    Cs_r=np.zeros((self.Nlevels),np.float64)
    z_r=np.zeros((len(sc_r),self.eta_rho,self.xi_rho),np.float64)

    h =self.depth

    hc=h.min() #np.min(h.min(),self.Tcline)
    
    for k in xrange(self.Nlevels):
        """
        Original vertical strectching function, Song and Haidvogel (1994).
        see: https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
        C(s) = (1 - b) * [SINH(s * a) / SINH(a)] +

             b * [-0.5 + 0.5 * TANH(a * (s + 0.5)) / TANH(0.5 * a)]
        """
        sc_r[k]=-1.0+(self.Nlevels-k-0.5)/self.Nlevels
        
        if self.vstretching==1:
            Ptheta=np.sinh(self.theta_s*sc_r[k])/np.sinh(self.theta_s)
            Rtheta=(np.tanh(self.theta_s*(sc_r[k]+0.5))/(2.0*np.tanh(0.5*self.theta_s)))-0.5

            if self.theta_s != 0.0:
                Cs_r[k]=(1-self.theta_b)*Ptheta+self.theta_b*Rtheta
            else:
                Cs_r[k]=sc_r[k]
                
            if k==0:
                Cs_r[k]=0

        if self.vstretching==2:
            """
            A. Shchepetkin new vertical stretching function.
            C(s) = [1.0 - COSH(theta_s * s)] / [COSH(theta_s) - 1.0]
            Adapted from set_scoord.F (ROMS source code)
            """
            Aweight=1.0
            Bweight=1.0

            if self.theta_s > 0.0:
                Csur=(1.0 - np.cosh(self.theta_s * sc_r[k]))/(np.cosh(self.theta_s) - 1.0)

                if self.theta_b > 0.0:

                    Cbot=np.sinh(self.theta_b)*(sc_r[k] + 1.0)/np.sinh(self.theta_b) - 1.0

                    Cweight=((sc_r[k]+1.0)**Aweight)*(1.0 + (Aweight/Bweight)*(1.0-(sc_r[k]+1.0)**Bweight))

                    Cs_r[k]=Cweight*Csur +(1.0-Cweight)*Cbot
                else:
                    Cs_r[k]=Csur
            else:
                Cs_r[k]=sc_r

    """
    TODO: FIXME hardcode variables should be read from file
    """
    zeta=None

    for k in xrange(len(sc_r)):
        if self.vstretching==1:
            z_r[k,:,:] =np.multiply(sc_r[k],hc) + np.subtract(h,hc)*Cs_r[k]
        
        if self.vstretching==2:
            z_r[k,:,:] =np.multiply(sc_r[k],hc) + np.subtract(h,hc)*Cs_r[k]

    self.z_r = z_r * self.mask_rho
    self.Cs_rho=Cs_r
    self.s_rho=sc_r
