# -*- coding: utf-8 -*-
"""
Reading csv files with ; for comments

This demonstrates how to read a csv file with comments.
I've assumed that comment lines start with a semicolon as
is the case with .nk files.
"""
import numpy as np
import matplotlib.pyplot as plt

wl = []
ndx = []
with open('Al.nk') as f:
    for line in f.readlines():
        if line[0]!=';':
            row=line.split();
            wl.append(float(row[0])/10);
            ndx.append(complex(row[1]+'+'+row[2]+'j'));
wl=np.array(wl)
ndx=np.array(ndx)
plt.semilogx(wl,ndx.real)
plt.title('Index of Refraction of Al')
plt.xlabel('wavelength, nm')
plt.ylabel('R(n)')
plt.figure()
plt.loglog(wl, ndx.imag)
plt.title('Index of Refraction of Al')
plt.xlabel('wavelength, nm')
plt.ylabel('I(n)')