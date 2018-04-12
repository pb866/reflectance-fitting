# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 01:20:57 2018

@author: rstur_000
"""
import refl
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
AlIndex = refl.Index('Al')
Al2O3Index = refl.Index('Al2O3')
npts=200
wl=np.linspace(50,200,npts)
ndx=AlIndex.at(wl)
ndx2=Al2O3Index.at(wl)
thetad=75
rfl1=[refl.matR(np.array([1,ndx2[i],ndx[i]]),np.array([0,3,0]),thetad,
                wl[i],0)[0] for i in range(npts)]
rfl2=[refl.matR(np.array([1,ndx[i]]),np.array([0,0]),thetad,
                 wl[i],0)[0] for i in range(npts)]
plt.figure()
plt.plot(wl,rfl1,wl,rfl2)
plt.legend(('Bare Al','3 nm Al2O3'))
plt.title('Aluminum Reflectance')
plt.xlabel('wavelength, nm')
plt.ylabel('fraction')
plt.show()