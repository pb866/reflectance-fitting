# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 08:00:09 2018

empirically computing uncertainty

@author: rturley
"""

from refl import Parratt, Index
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


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

def print_fit(label, popt, rms):
    print(label)
    print("  real(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[0],rms[0]))
    print("  imag(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[1],rms[1]))
    print("  thickness = {0:7.4f} +/- {1:.5f}".format(
            popt[2],rms[2]))
        
sigma = np.ones(np.size(thr))

points = 100

print("CONSTANT ERROR\n")
popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=False)
rms = np.sqrt(np.diag(pcov))
print_fit("Calculated uncertainty", popt, rms)
print_fit("Exact", p0, np.zeros(3))

sump = np.zeros(3)
sumsq = np.zeros(3)
for i in range(points):
    refn = (rfl*(1+sigmap*np.random.randn(np.size(thr)))+
        sigmac*np.random.randn(np.size(thr)))
    popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=False)
    sump = sump + popt
    sumsq = sumsq + popt**2
mean = sump/points
rms = np.sqrt(sumsq/points-mean**2)
print_fit("Actual RMS", mean, rms)