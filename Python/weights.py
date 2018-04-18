# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:00:13 2018

Variable weights on reflectance fits

@author: rturley
"""

from refl import Parratt, Index
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
rfl=[Parratt(n,t,thetad, wl) for thetad in thr]
# Add noise
sigmap = 0.05 # proportiional noise
sigmac = 1e-4 # constant noise
refn = (rfl*(1+sigmap*np.random.randn(np.size(thr)))+
        sigmac*np.random.randn(np.size(thr)))

p0 = np.array([alf3ndx.real, alf3ndx.imag, 18])
def f(thr, n, k, t):
    ndx = np.array([1+0j, n+k*1j, alndx, sio2ndx, sindx])
    th = np.array([0, t, 50, 1.6, 0])
    return [Parratt(ndx ,th, thetad, wl) for thetad in thr]

def print_fit(label, popt, pcov):
    print(label)
    print("  real(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[0],np.sqrt(pcov[0,0])))
    print("  imag(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[1],np.sqrt(pcov[1,1])))
    print("  thickness = {0:7.4f} +/- {1:.5f}".format(
            popt[2],np.sqrt(pcov[2,2])))
        
sigma = np.ones(np.size(thr))
popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=False)
print(pcov[0,0])
print_fit('exact', np.array([alf3ndx.real, alf3ndx.imag, 18]),
          np.zeros([3,3]))
print_fit('unweighted fits', popt, pcov)
ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sio2ndx, sindx])
thf = np.array([0, popt[2], 50, 1.6, 0])
rfit = [Parratt(ndf ,thf, thetad, wl) for thetad in thr]

plt.close('all')
plt.figure()
plt.semilogy(thr, refn, '.', thr, rfit, '-', thr, rfl, '--')
plt.legend(['data','fit','exact'])
plt.title('Unweighted Fit for AlF3/Al/SiO2/Si with Noise')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')

# Weighted fit
sigma=sigmap*np.array(rfl)
popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=True)
print_fit('prop weighted fits', popt, pcov)
ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sio2ndx, sindx])
thf = np.array([0, popt[2], 50, 1.6, 0])
rfit = [Parratt(ndf ,thf, thetad, wl) for thetad in thr]

plt.figure()
plt.semilogy(thr, refn, '.', thr, rfit, '-', thr, rfl, '--')
plt.legend(['data','fit','exact'])
plt.title('Prop Weighted Fit for AlF3/Al/SiO2/Si with Noise')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')
plt.show()

