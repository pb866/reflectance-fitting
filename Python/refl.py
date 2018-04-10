# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:52 2015

@author: rst
"""
import math
import numpy as np
import scipy.interpolate as intp

def matR(n, t, thetad, lam, sigma):
    """Reflectance and transmittance form a multilayer mirror
    
    Parameters
    ----------
    n : array of numbers
        index of refraction of the various layers, starting
        with vacuum.
    t : array of numbers
        thickness of the various layers, starting with vacuum.
        The thickness of the vacuum and substrate should be 0.
    thetad: number
        incident angle in degrees
    sigma: array of numbers
        roughness at each interace, starting with vacumm/top layer
    
    Returns
    -------
    tuple with reflectance and transmittance
    
    Example
    -------
    >>> matR(n, t, 45, 30.4, sigma)
    """
    ky=2*np.pi*np.cos(thetad*np.pi/180)/lam
    ki = lambda ni : 2*np.pi*ni/lam
    kzi = lambda ni : np.sqrt(ki(ni)**2-ky**2)
    Ci = lambda ni, di : np.exp(np.complex(0, kzi(ni))*di/2)
    def rmats(n1, n2, d1, d2):
        kz1=kzi(n1)
        kz2=kzi(n2)
        fs12=(kz1-kz2)/(kz1+kz2)
        gs12=2*kz1/(kz1+kz2)
        fs21=-fs12
        gs21=2*kz2/(kz1+kz2)
        C1=Ci(n1,d1)
        C2=Ci(n2,d2)
        return np.array([[gs21*C1*C2-fs21*fs12*C1*C2/gs12,
                          fs12*C1/(gs12*C2)],
                        [-fs21*C2/(gs12*C1),
                         1/(gs12*C1*C2)]])
    def rmatp(n1, n2, d1, d2):
        kz1=kzi(n1)
        kz2=kzi(n2)
        fp12=(n2**2*kz1-n1**2*kz2)/(n1**2*kz2+n2**2*kz1)
        gp12=2*n1*n2*kz1/(n1**2*kz2+n2**2*kz1)
        fp21=-fp12
        gp21=2*n1*n2*kz2/(n1**2*kz2+n2**2*kz1)
        C1=Ci(n1,d1)
        C2=Ci(n2,d2)
        return np.array([[gp21*C1*C2-fp21*fp12*C1*C2/gp12,
                          fp12*C1/(gp12*C2)],
                        [-fp12*C2/(gp12*C1),
                         1/(gp12*C1*C2)]])
    A=np.identity(2)
    B=A
    for i in range(n.size-1,1,-1):
        A=np.matmul(rmats(n[i], n[i-1], t[i], t[i-1]),A)
        B=np.matmul(rmatp(n[i], n[i-1], t[i], t[i-1]),B)
    ts=1/A[1,1]
    rs=ts*A[0,1]
    tp=1/B[1,1]
    rp=tp*B[0,1]
    percentS=fracs(lam)
    r=percentS*np.abs(rs)**2+(1-percentS)*np.abs(rp)**2
    # This assumes starting and ending materials are the same
    t=percentS*np.abs(ts)**2+(1-percentS)*np.abs(tp)**2
    return (r,t)

def Parratt(n1, n2, thetad):
    """Reflectance from a multilayer mirror
    
    Parameters
    ----------
    n1 : number
        index of the incident layer
    n2 : number
        index of the other layer
    thetad : number
        incident angle in degrees
        
    Returns
    -------
    ref: number
        reflectance of the layer for s polarization
        
    Example
    -------
    >>> Parratt(1,0.9+0.05j,15)
    0.015
    """
    theta=thetad*math.pi/180
    costt=(1-n1/n2*math.sin(theta))**0.5
    n1c=n1*math.cos(theta)
    n2c=n2*costt
    return abs((n1c-n2c)/(n1c+n2c))**2

ipts = np.array([[1.05685E+1,9.37888E-1],
                 [1.40446E+1,9.26536E-1],
                 [1.75270E+1,9.15094E-1],
                 [2.25761E+1,9.02016E-1],
                 [2.95528E+1,8.85600E-1],
                 [3.80743E+1,8.70842E-1],
                 [4.67841E+1,8.57696E-1],
                 [5.74863E+1,8.44550E-1],
                 [7.06822E+1,8.26362E-1],
                 [8.69072E+1,8.08173E-1],
                 [1.15690E+2,7.83378E-1],
                 [1.46820E+2,7.63555E-1],
                 [1.83342E+2,7.47070E-1],
                 [2.40104E+2,7.27294E-1],
                 [3.04907E+2,7.02429E-1],
                 [3.99133E+2,6.86013E-1],
                 [5.14002E+2,6.74616E-1],
                 [6.72413E+2,6.63242E-1],
                 [9.06564E+2,6.61998E-1],
                 [1.04406E+3,6.63885E-1]])
xfr=ipts[:,0]
yfr=ipts[:,1]
frfunc=intp.interp1d(np.log10(xfr),yfr,'cubic')
def fracs(lam):
    """fraction of s poplarization at the ALS
    
    Parameters
    ----------
    lam : number
        wavelength in nm
        
    Returns
    -------
    fr : number
        fraction of s polarization at this wavelength
        
    Example
    -------
    >>> fracs(25.1)
    0.9272135287221943
    """
    ev=np.log10(1239.8/lam)
    return (frfunc(ev)+1)/2