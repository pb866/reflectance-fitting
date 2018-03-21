# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 02:11:09 2015

@author: rstur_000
"""

# Example reading a CSV file from ALS
import csv
import matplotlib.pyplot as plt

def alsdata(filename):
    pos=[];
    detector=[];
    with open(filename, newline='') as f:
        header=f.readline(); # skip header
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            pos.append(float(row[0]));
            signal = float(row[1]);
            current = float(row[3]);
            detector.append(signal*current/500);
    return (pos, detector)

def savedata(filename, y, d):
    with open(filename,'w', newline='') as f:
        writer = csv.writer(f)
        for i in range(len(y)):
            writer.writerow([y[i],d[i]])

filename='y2o3001006.dat'
yp,det=alsdata(filename)
plt.plot(yp, det)
plt.xlabel('y position, mm')
plt.ylabel('detector signal')
# Next line will make it a log plot if included
# plt.yscale('log')
plt.title('Y Scan')
print('there are {0} data points'.format(len(yp)))
savedata('test.csv',yp,det)