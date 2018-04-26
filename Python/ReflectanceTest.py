# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 11:50:40 2018

@author: rturley
"""

import refl
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


log = refl.Log('../ALSdata','Feb2018')
runs = refl.Runs(log)

# dark current calculation
d7=np.mean(runs[101].diode)
d8=np.mean(runs[98].diode)
d9=np.mean(runs[99].diode)
d10=np.mean(runs[100].diode)
dark=(d7,d8,d9,d10)

# runs
s18a = refl.Reflectance((runs[119],),runs[97], dark)
plt.figure()
plt.semilogy(s18a.ang, s18a.rfl)
plt.title('Spectrum 18a')
s18b = refl.Reflectance((runs[120], runs[121]), runs[97], dark)
plt.figure()
plt.semilogy(s18b.ang, s18b.rfl)
plt.title('Spectrum 18b')

s18 = s18a.filter(2,18) + s18b.filter(18.05,80)
s18.plot()

# index objects
AlF3Index = refl.Index("AlF3")
AlIndex = refl.Index("Al")
Si3N4Index = refl.Index("Si3N4")
# interpolated values at 18 nm
alf3ndx = AlF3Index.at(s18.wavelength)
alndx = AlIndex.at(s18.wavelength)
si3n4ndx = Si3N4Index.at(s18.wavelength)

# fit function
def f(thr, n, k, tf, ta):
    ndx = np.array([1+0j, n+k*1j, alndx, si3n4ndx])
    th = np.array([0, tf, ta, 0])
    sigma = (1.0,0,0);
    return [refl.Parratt(ndx ,th, thetad, s18.wavelength, 0, sigma)
            for thetad in thr]
