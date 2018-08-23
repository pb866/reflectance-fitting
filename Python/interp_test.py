# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 10:31:11 2018

@author: rturley
"""

from pandas import read_table as pread
import numpy as np
import matplotlib.pyplot as plt

class Index:
    """ Index of refraction from volta
    
    Constructor Parameters
    ----------------------
    material : string
        material from base of .nk file on volta
    
    Method
    ------
    at(wavelength) : returns interpolated index at the given
        wavelength in nm. Wavelength can by an np.array
    
    Example
    -------
    alndx=Index('Al')
    lam=np.linspace(10,400,200)
    ndx=alndx.at(lam)
    """
    def __init__(self, material, lambda_min=0, lambda_max=100000):
        df=pread("http://volta.byu.edu/nk/"+material+".nk",
                           comment=';',delim_whitespace=True)
        val=df.values
        self.lam = [v for v in val[:,0] if (v>=lambda_min) and (v<=lambda_max)]
        self.n = [v[1] for v in val if (v[0]>=lambda_min) and (v[0]<=lambda_max)]
        self.k = [v[2] for v in val if (v[0]>=lambda_min) and (v[0]<=lambda_max)]
# =============================================================================
#         self.at=interp1d(lam,ndx,'cubic')
# =============================================================================

plt.close('all')
import scipy.interpolate as intr
ndx = Index('Al', lambda_min=350, lambda_max=500)
ndx2 = Index('Al2O3', lambda_min = 350, lambda_max = 500)
kintr = intr.PchipInterpolator(ndx2.lam, ndx2.n)
xx = np.linspace(ndx2.lam[1],ndx2.lam[-1], 300);
yy = kintr(xx);
sintr = intr.interp1d(ndx2.lam, ndx2.n, kind='cubic')
yyy = sintr(xx)
plt.figure()
plt.plot(ndx2.lam,ndx2.n,'.', xx, yy, '-', xx, yyy, '-')
plt.title('Al$_2$O$_3$ Index of Refraction')
plt.legend(('data','spline','PCHIP'))
plt.xlabel('wavelength, $\\AA$')
plt.ylabel('n')
plt.show()
plt.figure()
kintr = intr.PchipInterpolator(ndx.lam, ndx.k)
xx = np.linspace(355,495,300);
yy = kintr(xx);
sintr = intr.interp1d(ndx.lam, ndx.k, kind='cubic')
yyy = sintr(xx)
plt.semilogy(ndx.lam,ndx.k,'.',xx,yy,'-',xx,yyy,'-')
plt.title('Aluminum Index of Refraction')
plt.legend(('data','PCHIP','spline'))
plt.xlabel('wavelength, $\\AA$')
plt.ylabel('k')
plt.show()