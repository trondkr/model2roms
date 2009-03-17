
"""
Module with mathematically defined orthogonal grids
"""

# x, y = ROMS grid coordinates
# x = i, gives i-th rho-point (counting from 0 to L)
# similarly for j

from numpy import *
import numpy as np

# ------------------------------
# Useful constants
# ------------------------------
# sjekk griffies' bok for verdi
# 6367 => best fit to geopotential (Gill)
# 6371 => volume of sphere = vol. of earth (Gill)
# Sjekk hva Kate Hedstrom bruker
R_earth = 6378.137 # Earth radius   [km]  
pi      = np.pi
rad     = pi/180.0
deg     = 180.0/pi



# -----------------------------
# Useful functions
# -----------------------------



# Compute distances on a spherical earth with Haversine formula
def dist(lon0, lat0, lon1, lat1):
    """Computes distances on spherical earth, Haversine formula"""
    phi0    = lat0*rad
    phi1    = lat1*rad
    dphi    = phi1 - phi0
    dlambda = (lon1 - lon0) * rad
    a = sin(0.5*dphi)**2 + cos(phi0)*cos(phi1)*sin(0.5*dlambda)**2
    return 2 * R_earth * arctan2(sqrt(a), sqrt(1-a))


def ll2xyz(lon, lat):
    """Returns Cartesian coordinates x, y, z from longitude and latitude"""
    lambda_ = lon*rad
    cos_phi = cos(lat*rad)    
    x = cos(lambda_)*cos_phi
    y = sin(lambda_)*cos_phi
    z = sin(lat*rad)
    return x, y, z

def xyz2ll(x, y, z):
    """Returns longitude and latitude from cartesian coordinates,
       Projects onto unit sphere
    """
    r = sqrt(x*x + y*y + z*z)
    lat = np.arcsin(z/r)*deg
    lon = np.arctan2(y,x)*deg
    return lon, lat

# Missing function in ipython,
# arctan2 version is slightly faster
def arg(z):
    """Return arg s.t. z = abs(z)*exp(1j*arg(z))"""
    #return log(z/abs(z)).imag
    return arctan2(z.imag, z.real)


# Polar stereographic: sphere S^2 -> ekvator plane C

def pst(lon, lat):
    """Polar stereographic map onto complex ekvator plane"""
    return np.exp(1j*lon*rad)*np.tan(0.25*pi-0.5*lat*rad)

# Global inverse polarstereographic
# fra complex projective line to unit sphere
# Bug?
def pst_inv2(z1,z2=1):
    #z1 = complex(z1)
    #z2 = complex(z2)
    z1 = np.array(z1, dtype='complex128')
    z2 = np.array(z2, dtype='complex128')
    z1c = z1.conjugate()
    z2c = z2.conjugate()
    n1 = (z1*z1c).real
    n2 = (z2*z2c).real
    x = (z1*z2c).real / (n1+n2)
    y = (z1*z2c).imag / (n1+n2)
    z = (n2-n1) / (n1 + n2)
    return xyz2ll(x, y, z)

def pst_inv(z):
    """Inverse polar stereographic map"""
    r = abs(z)
    lat = 90.0 - 2*arctan(r)*deg
    lon = arg(z)*deg
    return lon, lat
        
def merc(lon, lat):
    """Mercator projection map"""
    return lon*rad + 1j*log(tan(0.25*pi+0.5*lat*rad))

def merc_inv(z):
    """Inverse mercator projection map"""
    lon = z.real*deg
    lat = 2*arctan(exp(z.imag))*deg - 90.0
    return lon, lat
    


# --- Square mercator grid ---

class square_mercator:
    def __init__(self, lon0, lat0, dlon):
        self.lon0 = lon0
        self.lat0 = lat0
        self.dlon = dlon

    def grid2ll(self, x, y):
        h = self.dlon*rad
        phi0 = self.lat0*rad
        return ( self.lon0 + x*self.dlon,
                 (0.5*pi - 2*arctan(tan(0.25*pi-0.5*phi0)*exp(-y*h)))*deg
               )

    def ll2grid(self, lon, lat):
        h       = self.dlon*rad
        lambda0 = self.lon0*rad
        phi0    = self.lat0*rad
        phi     = lat*rad
        return ( (lon*rad - lambda0)/h,
                  log( tan(0.25*pi-0.5*phi0) / tan(0.25*pi-0.5*phi) ) / h
               )

    def angle(self, x, y):
        return 0.0


# --- Polarstereografisk "feltfile" grid

class square_stereographic:
    def __init__(self, xp, yp, dx, ylon):
        self.xp   = xp         
        self.yp   = yp
        self.dx   = dx
        self.ylon = ylon

   
    def grid2ll(self, x, y):
        phi0    = 60.0*rad
        lambda0 = self.ylon*rad
        xp  = self.xp
        yp  = self.yp
        dx  = self.dx
        lon = lambda0 + arctan2(x - xp, yp - y)
        r   = dx * sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
        lat = 0.5*pi - 2*arctan(r / (R_earth*(1+sin(phi0))))
        return (lon*deg, lat*deg)
        
    def ll2grid(self, lon, lat):
        phi0    = 60.0*rad
        lambda0 = self.ylon*rad
        lambda_ = lon*rad
        phi     = lat*rad
        xp  = self.xp
        yp  = self.yp
        dx  = self.dx
        r = R_earth * cos(phi) * (1+sin(phi0)) / (1+sin(phi))
        x = xp + r*sin(lambda_ - lambda0)/dx
        y = yp - r*cos(lambda_ - lambda0)/dx
        return (x, y)

    def ll2grid2(self, lon, lat):
        phi0    = 60.0*rad
        lambda0 = self.ylon*rad
        #lambda_ = lon*rad
        #phi     = lat*rad
        xp  = self.xp
        yp  = self.yp
        dx  = self.dx
        scale = R_earth*(1+sin(phi0))/dx
        w = (xp + 1j*yp) - scale*exp(-1j*lambda0)*1j*pst(lon,lat)
        return w.real, w.imag

    def angle(self, x, y):
        return arctan2(x-self.xp, y-self.yp)


# Polatsterographic square polar coordinates
# Bedre navn ?
# Gitt ved data:
# lengde, bredde for sentrum  S
# lengde, bredde for origo    O
# vinkelsteglengde, dangle    

class stereographic_wedge:
    def __init__(self, lonS, latS, lonO, latO, dangle):
        self.lonS = lonS
        self.latS = latS
        self.lonO = lonO
        self.latO = latO
        self.dangle = dangle

    def grid2ll(self, x, y):
        # Forenkler formlene med kompleks aritmetikk
        phi0    = 60.0*rad
        lambdaS = self.lonS*rad
        phiS    = self.latS*rad
        lambdaO = self.lonO*rad
        phiO    = self.latO*rad
        alpha   = self.dangle*rad

        rS = R_earth * cos(phiS) * (1+sin(phi0)) / (1+sin(phiS))
        rO = R_earth * cos(phiO) * (1+sin(phi0)) / (1+sin(phiO))
        j = complex(0, 1)
        zS = rS*exp(j*lambdaS)
        zO = rO*exp(j*lambdaO)
        d0 = abs(zS+zO)
                
        theta = y * alpha
        d = d0 * exp(alpha*x)

        #print "d0, d = ", d0, d
        #print "theta = ", theta

        zP = zS + (zO-zS)*d*exp(j*theta)/d0
        rP = abs(zP)
        lambdaP = arctan2(zP.imag, zP.real)
        phiP = 0.5*pi - 2*arctan(rP / (R_earth*(1+sin(phi0))))
        return (lambdaP*deg, phiP*deg)

    def ll2grid(self, lon, lat):
        lambda_ = lon*rad
        phi     = lat*rad
        phi0    = 60.0*rad
        lambdaS = self.lonS*rad
        phiS    = self.latS*rad
        lambdaO = self.lonO*rad
        phiO    = self.latO*rad
        alpha   = self.dangle*rad
        
        r =  R_earth * cos(phi) * (1+sin(phi0)) / (1+sin(phi))
        rS = R_earth * cos(phiS) * (1+sin(phi0)) / (1+sin(phiS))
        rO = R_earth * cos(phiO) * (1+sin(phi0)) / (1+sin(phiO))
        j = complex(0, 1)
        w  = r  * exp(j*lambda_)
        wS = rS * exp(j*lambdaS)
        wO = rO * exp(j*lambdaO)

        # exp(a(x+j*y)) = (w-wS)/(wO-wS)
        z1 = (w - wS)/(wO - wS)
        #x1 = z1.real
        #y1 = z1.imag
        #x  = log(sqrt(x1*x1 + y1*y1))/alpha
        #y  = arctan2(y1,x1)/alpha
        z = log(z1)/alpha
       
        return (z.real, z.imag)

# ------

# utility functions for transformed mercator




class transformed_mercator:

    # Note that slim and nlim are not used for mercator   
    def __init__(self, lon_a, lat_a, lon_b, lat_b, wlim, elim, jres, ires):
        self.a = pst(lon_a, lat_a)
        self.b = pst(lon_b, lat_b)

        xa, ya, za = ll2xyz(lon_a, lat_a)
        xb, yb, zb = ll2xyz(lon_b, lat_b)
        xc = xa + xb
        yc = ya + yb
        zc = za + zb
        self.c = pst(*xyz2ll(xc, yc, zc))

        self.wlim = wlim
        self.elim = elim
        self.jres = jres
        self.dm = (elim-wlim)/float(jres)
        self.slim = -0.5*ires*self.dm
        self.nlim = 0.5*ires*self.dm

    def grid2ll(self, x, y):
        # Identical to Paul's use of Bentsen fortran code
        a = self.a
        b = self.b
        c = self.c
        
        lon1 = self.elim - (y-0.5)*self.dm
        lat1 = - 90.0 + 2*arctan(exp((self.slim+(x-0.5)*self.dm)*rad))*deg
        w = pst(lon1, lat1)
        z = (-b*w*(c-a)+a*(c-b))/(-w*(c-a)+(c-b))
        return pst_inv(z)

    # OK, er invers
    def ll2grid(self, lon, lat):
        
        a = self.a
        b = self.b
        c = self.c
        
        z = pst(lon, lat)
        w = (z-a)*(c-b)/((z-b)*(c-a))
        lon1, lat1 = pst_inv(w)
        x = 0.5 + (-self.slim + log(tan(0.25*pi + 0.5*lat1*rad))*deg)/self.dm
        y = 0.5 + (self.elim - lon1)/self.dm

        return x, y

    def ll2grid2(self, lon, lat):
        # OK denn er riktig
        a = self.a
        b = self.b
        c = self.c
        
        z = pst(lon, lat)
        w = (z-a)*(c-b)/((z-b)*(c-a))
        lon1, lat1 = pst_inv(w)
        wm = 0.5*(1+1j) -     \
           1j*(merc(lon1,lat1)*deg - complex(self.elim, self.slim))/self.dm
        return wm.real, wm.imag



        
# ---------------------------------------

# lon0 = minimum longitude
# lat0 = minimum latitiude
# dlon = step-lengde i longitude
# dlat = step-lengde i latitude

class spherical:
    def __init__(self, lon0, lat0, dlon, dlat):
        self.lon0 = lon0
        self.lat0 = lat0
        self.dlon = dlon
        self.dlat = dlat

    def grid2ll(self, x, y):
        return (self.lon0 + x*self.dlon, self.lat0 + y*self.dlat)

    def ll2grid(self, lon, lat):
        return ((lon - self.lon0)/self.dlon, (lat-self.lat0)/self.dlat)
    
        
    
    
