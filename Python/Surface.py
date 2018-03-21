# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 10:37:00 2015

@author: Steve Turley
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

    
class Surface:
    """class for representing a 1d surface
    
    Parameters
    ----------
    p2: string
        name of the surface
    
    Constructor Parameters
    ----------------------
    n1 : number
        index of the incident layer
    n2 : number
        index of the other layer
    thetad : number
        incident angle in degrees
                
    Example
    -------
    >>> refl(1,0.9+0.05j,15)
    0.015
    """
    
    def __init__(self, sigma, length, npts):
        """Constructor to initialize surface
        
        Parameters
        ----------
        sigma : number
            rms surface roughness
        length : number
            length of surface
        npts : integer
            number of points in surface
        
        Returns
        -------
        nothing
        
        Example
        -------
        >>> surf = Surface(0.1,5.0,100)
            creates surface
        """
        self.xp=np.linspace(0,length,npts)
        self.yp=stats.norm.rvs(0,sigma,npts)
    
    def plot(self):
        plt.plot(self.xp,self.yp)
