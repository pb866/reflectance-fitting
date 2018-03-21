# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:52 2015

@author: rst
"""
import math
def refl(n1, n2, thetad):
    """Reflectance for s polarization
    
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
    res: number
        reflectance of the layer for s polarization
        
    Example
    -------
    >>> refl(1,0.9+0.05j,15)
    0.015
    """
    theta=thetad*math.pi/180
    costt=(1-n1/n2*math.sin(theta))**0.5
    n1c=n1*math.cos(theta)
    n2c=n2*costt
    return abs((n1c-n2c)/(n1c+n2c))**2