# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:52 2015

@author: rst
"""
import math
import numpy as np
import scipy.interpolate as intp

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