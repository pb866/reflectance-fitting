# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 11:50:40 2018

@author: rturley
"""

import refl
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# log = refl.Log('../ALSdata','Feb2018') # NSF computer
log = refl.Log() # BYU Computer
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
    sigma = np.array([1.0,0,0]);
    return [refl.Parratt(ndx ,th, thetad, s18.wavelength, 0, sigma)
            for thetad in thr]

alpha = 0.2
beta = 0.097
gamma = 5e-5
sigmap = np.array([alpha*np.exp(-beta*th) for th in s18.ang])
sigmac = gamma
sigmaw = np.sqrt(sigmap**2+sigmac**2)
p0=np.array([0.95, 0.025, 8, 26])
popt, pcov = curve_fit(f, s18.ang, s18.rfl, p0, sigmaw,
                       absolute_sigma=False)
print('n = '+str(round(popt[0],4))+"+/-"+
      str(round(np.sqrt(pcov[0,0]),4)))
print('k = '+str(round(popt[1],4))+"+/-"+
      str(round(np.sqrt(pcov[1,1]),2)))
print('tf = '+str(round(popt[2],2))+"+/-"+
      str(round(np.sqrt(pcov[2,2]),2)))
print('ta = '+str(round(popt[3],4))+"+/-"+
      str(round(np.sqrt(pcov[3,3]),2)))
rfit = f(s18.ang, popt[0], popt[1], popt[2], popt[3])
plt.figure()
plt.semilogy(s18.ang, s18.rfl,'.',s18.ang, rfit, '-',s18.ang, sigmaw,'-r')
plt.title('Data and Fit')
plt.xlabel('grazing angle, deg')
plt.ylabel('reflectance')
plt.legend(['data', 'fit','sigma'])

# refit with weights proportional to the data
sigmap = alpha*np.abs(s18.rfl)
sigmaw = np.sqrt(sigmap**2+sigmac**2)
popt, pcov = curve_fit(f, s18.ang, s18.rfl, p0, sigmaw,
                       absolute_sigma=False)
print('n = '+str(round(popt[0],4))+"+/-"+
      str(round(np.sqrt(pcov[0,0]),4)))
print('k = '+str(round(popt[1],4))+"+/-"+
      str(round(np.sqrt(pcov[1,1]),2)))
print('tf = '+str(round(popt[2],2))+"+/-"+
      str(round(np.sqrt(pcov[2,2]),2)))
print('ta = '+str(round(popt[3],4))+"+/-"+
      str(round(np.sqrt(pcov[3,3]),2)))
rfit = f(s18.ang, popt[0], popt[1], popt[2], popt[3])
plt.figure()
plt.semilogy(s18.ang, s18.rfl,'.',s18.ang, rfit, '-')
plt.title('Data and Fit')
plt.xlabel('grazing angle, deg')
plt.ylabel('reflectance')
plt.legend(['data', 'fit'])

# Check \chi^2
chi2 = np.sum(((s18.rfl-rfit)/sigmaw)**2)
print('chi**2 = '+str(round(chi2,1))+" npts = "+str(np.size(sigmaw)))

# Add Al2O3 layer
Al2O3Index = refl.Index("Al2O3")
al2o3ndx = Al2O3Index.at(s18.wavelength)
def f2(thr, n, k, tf, ta, to):
    ndx = np.array([1+0j, n+k*1j, al2o3ndx, alndx, si3n4ndx])
    th = np.array([0, tf, to, ta, 0])
    sigma = np.array([1.0,0,0,0]);
    return [refl.Parratt(ndx ,th, thetad, s18.wavelength, 0, sigma)
            for thetad in thr]
p0=np.array([0.95, 0.025, 8, 26, 1])
popt, pcov = curve_fit(f2, s18.ang, s18.rfl, p0, sigmaw,
                       absolute_sigma=False)
print('n = '+str(round(popt[0],4))+"+/-"+
      str(round(np.sqrt(pcov[0,0]),4)))
print('k = '+str(round(popt[1],4))+"+/-"+
      str(round(np.sqrt(pcov[1,1]),2)))
print('tf = '+str(round(popt[2],2))+"+/-"+
      str(round(np.sqrt(pcov[2,2]),2)))
print('ta = '+str(round(popt[3],4))+"+/-"+
      str(round(np.sqrt(pcov[3,3]),2)))
print('to = '+str(round(popt[4],4))+"+/-"+
      str(round(np.sqrt(pcov[4,4]),4)))
rfit = f2(s18.ang, popt[0], popt[1], popt[2], popt[3], popt[4])
plt.figure()
plt.semilogy(s18.ang, s18.rfl,'.',s18.ang, rfit, '-')
plt.title('Data and Fit with Oxide Layer')
plt.xlabel('grazing angle, deg')
plt.ylabel('reflectance')
plt.legend(['data', 'fit','sigma'])

# Fit Al index of refraction
def f3(thr, nf, kf, na, ka, tf, ta):
    ndx = np.array([1+0j, nf+kf*1j, na+ka*1j, si3n4ndx])
    th = np.array([0, tf, ta, 0])
    sigma = np.array([1.0,0,0]);
    return [refl.Parratt(ndx ,th, thetad, s18.wavelength, 0, sigma)
            for thetad in thr]
p0=np.array([0.95, 0.025, alndx.real, alndx.imag, 8, 26])
popt, pcov = curve_fit(f3, s18.ang, s18.rfl, p0, sigmaw,
                       absolute_sigma=False)
print('nf = '+str(round(popt[0],4))+"+/-"+
      str(round(np.sqrt(pcov[0,0]),4)))
print('kf = '+str(round(popt[1],4))+"+/-"+
      str(round(np.sqrt(pcov[1,1]),4)))
print('na = '+str(round(popt[2],4))+"+/-"+
      str(round(np.sqrt(pcov[2,2]),4)))
print('ka = '+str(round(popt[3],4))+"+/-"+
      str(round(np.sqrt(pcov[3,3]),4)))
print('tf = '+str(round(popt[4],4))+"+/-"+
      str(round(np.sqrt(pcov[4,4]),4)))
print('ta = '+str(round(popt[5],4))+"+/-"+
      str(round(np.sqrt(pcov[5,5]),4)))
rfit = f3(s18.ang, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
plt.figure()
plt.semilogy(s18.ang, s18.rfl,'.',s18.ang, rfit, '-')
plt.title('Data and Fit Varyin Al Index')
plt.xlabel('grazing angle, deg')
plt.ylabel('reflectance')
plt.legend(['data', 'fit','sigma'])
