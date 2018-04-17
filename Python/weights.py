# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:00:13 2018

Variable weights on reflectance fits

@author: rturley
"""

from refl import matR, Index
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

AlIndex = Index("Al")
SiIndex = Index("Si")
SiO2Index = Index("SiO2")
AlF3Index = Index("AlF3")

wl = 12
alndx = AlIndex.at(wl)
sindx = SiIndex.at(wl)
sio2ndx = SiO2Index.at(wl)
alf3ndx = AlF3Index.at(wl)

n=np.array([1+0j, alf3ndx, alndx, sio2ndx, sindx])
t=np.array([0, 18, 50, 1.6, 0])
thr = np.linspace(0.5,80,160)
refl=[matR(n,t,thetad, wl, 0.0)[0] for thetad in thr]
# Add noise
refn = refl*(1+0.05*np.random.randn(np.size(thr)))+1e-4*np.random.randn(np.size(thr))

p0 = np.array([alf3ndx.real, alf3ndx.imag, 18])
def f(thr, n, k, t):
    ndx = np.array([1+0j, n+k*1j, alndx, sio2ndx, sindx])
    th = np.array([0, t, 50, 1.6, 0])
    return [matR(ndx ,th, thetad, wl, 0.0)[0] for thetad in thr]

sigma = np.ones(np.size(thr))
popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=False)
print(popt)
ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sio2ndx, sindx])
thf = np.array([0, popt[2], 50, 1.6, 0])
rfit = [matR(ndf ,thf, thetad, wl, 0.0)[0] for thetad in thr]

plt.close('all')
plt.figure()
plt.semilogy(thr, refn, '.', thr, rfit, '-', thr, refl, '--')
plt.legend(['data','fit','exact'])
plt.title('Fit for AlF3/Al/SiO2/Si with Noise')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')
plt.show()

