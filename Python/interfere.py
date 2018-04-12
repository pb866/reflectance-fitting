# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:03:26 2018

@author: rturley
"""

import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,3,200)
y1=np.sin(2*np.pi*x)
y2=0.9*np.sin(2*np.pi*x)
y3=-0.9*y1
plt.close('all')
plt.figure()
plt.plot(x,y1,x,y2)
plt.title('Constructive Interference')
plt.figure()
plt.plot(x,y1,x,y3)
plt.title('Destructive Interference')
plt.show()