# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:52 2015

@author: rst
"""
import math
import numpy as np
import scipy.interpolate as intp

def Parratt(n, x, thetad, lam, sigma):
    """Reflectance from a multilayer mirror
    
    Parameters
    ----------
    n : np.array of index of refractions for stack starting with 
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
# =============================================================================
# function R=Parratt(n,x,theta,fractions,lambda,sigma)
# % calculate the reflectance a multilayer stack
# % Written by Steve Turley, 6 July 2006
# fractionp=1-fractions; % fraction of p-polarized light
# S=sqrt(n.^2-cosd(theta)^2);
# k=2*pi/lambda;
# C=exp(i*2*S.*x*k);
# rs=0;
# rp=0;
# 
# qz=k*sind(theta); %Debye-Waller rougness correction
# eta=exp(-2*qz^2*sigma.^2); %Debye-Waller rougness correction
# 
# %qz=k*S; %Nove-Croce rougness correction
# for m=length(n)-1:-1:1;
#     %eta=exp(-2*qz(m)*qz(m+1)*sigma(m)^2); %Nove-Croce rougness correction
#     fs = (S(m)-S(m+1))/(S(m)+S(m+1));
#     fp = (n(m+1)^2*S(m)-n(m)^2*S(m+1))/(n(m+1)^2*S(m)+n(m)^2*S(m+1));
#     rs=C(m)*(fs*eta(m)+rs*eta(m)^2)/(1+fs*rs*eta(m));
#     rp=C(m)*(fp*eta(m)+rp*eta(m)^2)/(1+fp*rp*eta(m));
# end
# R=(fractionp*abs(rp)^2+fractions*abs(rs)^2);
# return;
# =============================================================================
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