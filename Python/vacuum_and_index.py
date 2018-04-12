# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 14:22:13 2018

@author: rturley
"""
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
plt.figure()
lam=np.linspace(0,3,200)
vac=np.cos(2*np.pi*lam)
npx=np.exp(-0.7*lam)*np.cos(1.5*np.pi*lam)
plt.plot(lam,vac,lam,npx)
plt.legend(('vacuum','with complex index'))
plt.title('Effect of Index of Refraction')
plt.show()