# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 05:52:43 2018

@author: rstev_000
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:00:13 2018

Variable weights on reflectance fits. Dealing with systematic
errors

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
    ndx = np.array([1+0j, n+k*1j, alndx, sindx])
    th = np.array([0, t, 50, 0])
    return [Parratt(ndx ,th, thetad, wl) for thetad in thr]

def print_fit(label, popt, pcov):
    print(label)
    print("  real(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[0],np.sqrt(pcov[0,0])))
    print("  imag(n)   = {0:7.4f} +/- {1:.5f}".format(
            popt[1],np.sqrt(pcov[1,1])))
    print("  thickness = {0:7.4f} +/- {1:.5f}".format(
            popt[2],np.sqrt(pcov[2,2])))
        
ndm = np.array([1, alf3ndx, alndx, sindx])
thm = np.array([0, 18, 50, 0])
refm = [Parratt(ndm, thm, thetad, wl) for thetad in thr]

sigma = np.ones(np.size(thr))
popt, pcov = curve_fit(f, thr, rfl, p0, sigma, absolute_sigma=False)
print_fit('exact', np.array([alf3ndx.real, alf3ndx.imag, 18]),
          np.zeros([3,3]))
print_fit('unweighted fits', popt, pcov)
ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sindx])
thf = np.array([0, popt[2], 50, 0])
rfit = [Parratt(ndf ,thf, thetad, wl) for thetad in thr]

plt.close('all')
plt.figure()
plt.semilogy(thr, rfl, '-', thr, refm, '-')
plt.title('Reflectance With and Without SiO2')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')
plt.legend(['with SiO2', 'without SiO2'])

plt.figure()
plt.semilogy(thr, rfl, '.', thr, rfit, '-')
plt.title('Fit Using Constant Weights')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')
plt.legend(['data','fit'])

res = 30000*(np.array(rfit)-np.array(rfl))/sigma
plt.figure()
plt.plot(thr, res, '.')
plt.title('Unweighted Fit with Systematic Error')
plt.xlabel('angle, degrees')
plt.ylabel('weighted residual')

# Proportionally Weighted fit
sigma=sigmap*np.array(rfl)
popt, pcov = curve_fit(f, thr, rfl, p0, sigma, absolute_sigma=False)
print_fit('prop weighted fits', popt, pcov)
ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sindx])
thf = np.array([0, popt[2], 50, 0])
rfit = [Parratt(ndf ,thf, thetad, wl) for thetad in thr]

plt.figure()
plt.semilogy(thr, rfl, '.', thr, rfit, '-')
plt.legend(['data','fit'])
plt.title('Prop Weighted Fit for Bad Model')
plt.xlabel('angle, degrees')
plt.ylabel('reflectance')

res = (np.array(rfit)-np.array(rfl))/sigma
plt.figure()
plt.plot(thr, res, '.')
<<<<<<< HEAD
plt.title('Normalized Residual for Proportionally Weighted Fit')
plt.xlabel('angle, degrees')
plt.ylabel('residual')

=======
plt.title('Proportionally Weighted Fit with Systematic Error')
plt.xlabel('angle, degrees')
plt.ylabel('weighted residual')

plt.show()
#
#
>>>>>>> 3a35739dadaae29ec46abf997a80ff30d379e9cf
## Combined Weighted fit
#sigma = np.sqrt((sigmap*np.array(rfl))**2+sigmac**2)
#popt, pcov = curve_fit(f, thr, refn, p0, sigma, absolute_sigma=True)
#print_fit('combined weighted fits', popt, pcov)
#ndf = np.array([1+0j, popt[0]+popt[1]*1j, alndx, sio2ndx, sindx])
#thf = np.array([0, popt[2], 50, 1.6, 0])
#rfit = [Parratt(ndf ,thf, thetad, wl) for thetad in thr]
#
#plt.figure()
#plt.semilogy(thr, refn, '.', thr, rfit, '-', thr, rfl, '--')
#plt.legend(['data','fit','exact'])
#plt.title('Combined Weighted Fit for AlF3/Al/SiO2/Si with Noise')
#plt.xlabel('angle, degrees')
#plt.ylabel('reflectance')
#
#res = (rfit-refn)/sigma
#plt.figure()
#plt.plot(thr, res, '.')
#plt.title('Residual for Combined Weighted Fit')
#plt.xlabel('angle, degrees')
#plt.ylabel('residual')
#
plt.show()
