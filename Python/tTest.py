# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 12:00:20 2018

@author: rstur_000
"""
import numpy as np
# import matplotlib.pyplot as plt

noise=[0.712705059,
    1.037487891,
    0.292159257,
    -0.71952627,
    -1.187286216]
m=15
b=3
x=np.linspace(1,5,5)
y=noise+m*x+b

p, cov = np.polyfit(x,y,1,full=False, cov=True)
pfit = np.poly1d(p)
xp = np.linspace(1,5,100)
# plt.plot(x,y,'.',xp,pfit(xp),'-')
std = np.sqrt(np.diag(cov))
print("slope     = {:8.5f} +/- {:.5f}".format(p[0],std[0]))
print("intercept = {:8.5f} +/- {:.5f}".format(p[1],std[1]))

# Try something else
p, S, mu = np.polyfit(x, y, 1)
y, delta = np.polyval(p, x, S, mu)
print("y=",y)
print("delta=",delta)

# Test these variances, they disagree with everyone else
sump = np.zeros(2)
sumsq = np.zeros(2)
sumstd = np.zeros(2)
trials=100
for i in range(100):
    noise = np.random.normal(scale=0.5,size=5)
    y = noise + m*x + b
    p, cov = np.polyfit(x,y,1,full=False, cov=True)
    sump += p
    sumsq += p**2
    sig = np.sqrt(np.diag(cov))
    sumstd += sig
meanp = sump/trials
meansq = sumsq/trials
meanstd = sumstd/trials
calcstd = np.sqrt(meansq-meanp**2)
print("calculted std = {:.5f}, {:.5f}".format(calcstd[0],calcstd[1]))
print("mean std: ",meanstd)
print("meanstd/calcstd = ", meanstd/calcstd)
