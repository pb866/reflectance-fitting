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
ndx = Index('Al', lambda_min=350, lambda_max=500)
# ndx2 = Index('Al2O3', lambda_min = 10, lambda_max = 1000)
# =============================================================================
# npts=200
# wl=np.linspace(10,1000,npts)
# ndx=AlIndex.at(wl)
# ndx2=Al2O3Index.at(wl)
# thetad=75
# =============================================================================
# =============================================================================
# plt.figure()
# plt.semilogy(ndx.lam, ndx.n,'.',ndx2.lam,ndx2.n,'.')
# plt.legend(('Al','Al2O3'))
# plt.title('Index of Refraction')
# plt.xlabel('wavelength, nm')
# plt.ylabel('n')
# plt.show()
# plt.figure()
# =============================================================================
import scipy.interpolate as intr
kintr = intr.PchipInterpolator(ndx.lam, ndx.k)
xx = np.linspace(355,495,300);
yy = kintr(xx);
sintr = intr.interp1d(ndx.lam, ndx.k, kind='cubic')
yyy = sintr(xx)
plt.semilogy(ndx.lam,ndx.k,'.',xx,yy,'-',xx,yyy,'-')
plt.title('Aluminum Index of Refraction')
plt.legend(('data','PCHIP','spline'))
plt.xlabel('wavelength, nm')
plt.ylabel('k')
plt.show()