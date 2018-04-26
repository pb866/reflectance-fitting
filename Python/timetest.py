# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 07:43:25 2018

@author: rstur_000
"""
import time
import numpy as np
el=100000
st=time.time()
l=[]
for i in range(el):
    l.append(i)
nl=np.array(l)
en=time.time()
print("array conversion took "+str(en-st)+" seconds")

st=time.time();
l=np.array([])
for i in range(el):
    l=np.append(l,i)
en=time.time()
print("np append took "+str(en-st)+" seconds")