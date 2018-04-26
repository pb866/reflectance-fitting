# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:45:52 2015

@author: rst
"""
#
# As you're developing this, you'll need to reimport this as it is
# updated. The following code fragment will do that for you:
# import importlib
# importlib.reload(refl)
#
# Testing: matR checked against Matlab refl.m code. Answers a sub-answers
# agreed to within 6 significant digits. I attribute the difference to
# interpolation differences.

import numpy as np
import scipy.interpolate as intp

def matR(n, t, thetad, lam, sigma=0):
    """Reflectance and transmittance form a multilayer mirror
    
    Parameters
    ----------
    n : array of numbers
        index of refraction of the various layers, starting
        with vacuum.
    t : array of numbers
        thickness of the various layers, starting with vacuum.
        The thickness of the vacuum and substrate should be 0.
    thetad: number
        incident angle in degrees
    sigma: array of numbers
        roughness at each interace, starting with vacumm/top layer
    
    Returns
    -------
    tuple with reflectance and transmittance
    
    Example
    -------
    >>> matR(n, t, 45, 30.4, sigma)
    """
    ky=2*np.pi*np.cos(thetad*np.pi/180)/lam
    ki = lambda ni : 2*np.pi*ni/lam
    kzi = lambda ni : np.sqrt(ki(ni)**2-ky**2)
    Ci = lambda ni, di : np.exp(np.complex(0, kzi(ni))*di/2)
    def rmats(n1, n2, d1, d2):
        kz1=kzi(n1)
        kz2=kzi(n2)
        fs12=(kz1-kz2)/(kz1+kz2)
        gs12=2*kz1/(kz1+kz2)
        fs21=-fs12
        gs21=2*kz2/(kz1+kz2)
        C1=Ci(n1,d1)
        C2=Ci(n2,d2)
        return np.array([[gs21*C1*C2-fs21*fs12*C1*C2/gs12,
                          fs12*C1/(gs12*C2)],
                        [-fs21*C2/(gs12*C1),
                         1/(gs12*C1*C2)]])
    def rmatp(n1, n2, d1, d2):
        kz1=kzi(n1)
        kz2=kzi(n2)
        fp12=(n2**2*kz1-n1**2*kz2)/(n1**2*kz2+n2**2*kz1)
        gp12=2*n1*n2*kz1/(n1**2*kz2+n2**2*kz1)
        fp21=-fp12
        gp21=2*n1*n2*kz2/(n1**2*kz2+n2**2*kz1)
        C1=Ci(n1,d1)
        C2=Ci(n2,d2)
        return np.array([[gp21*C1*C2-fp21*fp12*C1*C2/gp12,
                          fp12*C1/(gp12*C2)],
                        [-fp21*C2/(gp12*C1),
                         1/(gp12*C1*C2)]])
    A=np.identity(2)
    B=A
    for i in range(n.size-1,0,-1):
        A=np.matmul(rmats(n[i-1], n[i], t[i-1], t[i]),A)
        B=np.matmul(rmatp(n[i-1], n[i], t[i-1], t[i]),B)
#        print("When i="+str(i))
#        print("   A="+str(A))
#        print("   B="+str(B))
    ts=1/A[1,1]
    rs=ts*A[0,1]
    tp=1/B[1,1]
    rp=tp*B[0,1]
    try:
        percentS=fracs(lam)
    except ValueError:
        percentS=1.0
    r=percentS*np.abs(rs)**2+(1-percentS)*np.abs(rp)**2
    # This assumes starting and ending materials are the same
    t=percentS*np.abs(ts)**2+(1-percentS)*np.abs(tp)**2
    return (r,t)

def Parratt(n, x, thetad, lam, fractions=0, sigma=0):
    """Reflectance from a multilayer mirror
    
    Parameters
    ----------
    n : np.array of index of refractions for stack starting with 
        index of the incident layer
    x : np.array of layer thicknesses for stack starting with
        the thickness of the incident layer (usually vacuum).
        The vacuum and substrate thicknesses should be set to
        0.
    thetad : number
        incident angle in degrees
    lam : number
        wavelength in nanometers
    fractions : number
        fraction of light with s polarization. If zero, this is calculated
        using Gullikson formula for synchrotron
    sigma : array of number
        interface roughness in nm
        
    Returns
    -------
    ref: number
        reflectance of the layer for s polarization
        
    Example
    -------
    >>> lam = 15
    >>> AlIndex = refl.Index('Al')
    >>> alndx = AlIndex.at(lam)
    >>> SiO2Index = refl.index('SiO2')
    >>> sio2ndx = SiO2Index.at(lam)
    >>> n = np.array([1, alndx, sio2ndx])
    >>> t = np.array([0, 20, 0])
    >>> thetad = 20
    >>> Parratt(n, t, thetad, lam)
    0.012606490280013215
    """
    if fractions==0:
        fractions = fracs(lam)
    fractionp = 1-fractions
    S = np.sqrt(n**2-np.cos(thetad*np.pi/180)**2)
    k = 2*np.pi/lam
    C = np.exp(2j*S*x*k)
    rs = 0
    rp = 0
    
    qz = k*np.sin(thetad*np.pi/180) # for Debye-Waller correction
    if sigma == 0:
        sigma = np.zeros(n.size-1)
    eta = np.exp(-2*qz**2*sigma**2) # Debye-Waller roughness correction
    
    for m in range(n.size-1, 0,-1):
        fs = (S[m-1]-S[m])/(S[m-1]+S[m])
        fp = ((n[m]**2 * S[m-1] - n[m-1]**2 * S[m])/
              (n[m]**2 * S[m-1] + n[m-1]**2 * S[m]))
        rs = C[m-1]*(fs*eta[m-1]+rs*eta[m-1]**2)/(1+fs*rs*eta[m-1])
        rp = C[m-1]*(fp*eta[m-1]+rp*eta[m-1]**2)/(1+fp*rp*eta[m-1])
    return fractionp*np.abs(rp)**2+fractions*np.abs(rs)**2
     
# =============================================================================
# function R=Parratt(n,x,theta,fractions,lambda,sigma)
# % calculate the reflectance a multilayer stack
# % Written by Steve Turley, 6 July 2006
# fractionp=1-fractions; % fraction of p-polarized light
# S=sqrt(n.^2-cosd(theta)^2);
# k=2*pi/lambda;
# C=exp(i*2*S.*x*k);
# rs=0;
# rp=0;
# 
# qz=k*sind(theta); %Debye-Waller rougness correction
# eta=exp(-2*qz^2*sigma.^2); %Debye-Waller rougness correction
# 
# %qz=k*S; %Nove-Croce rougness correction
# for m=length(n)-1:-1:1;
#     %eta=exp(-2*qz(m)*qz(m+1)*sigma(m)^2); %Nove-Croce rougness correction
#     fs = (S(m)-S(m+1))/(S(m)+S(m+1));
#     fp = (n(m+1)^2*S(m)-n(m)^2*S(m+1))/(n(m+1)^2*S(m)+n(m)^2*S(m+1));
#     rs=C(m)*(fs*eta(m)+rs*eta(m)^2)/(1+fs*rs*eta(m));
#     rp=C(m)*(fp*eta(m)+rp*eta(m)^2)/(1+fp*rp*eta(m));
# end
# R=(fractionp*abs(rp)^2+fractions*abs(rs)^2);
# return;
# =============================================================================

ipts = np.array([[1.05685E+1,9.37888E-1],
                 [1.40446E+1,9.26536E-1],
                 [1.75270E+1,9.15094E-1],
                 [2.25761E+1,9.02016E-1],
                 [2.95528E+1,8.85600E-1],
                 [3.80743E+1,8.70842E-1],
                 [4.67841E+1,8.57696E-1],
                 [5.74863E+1,8.44550E-1],
                 [7.06822E+1,8.26362E-1],
                 [8.69072E+1,8.08173E-1],
                 [1.15690E+2,7.83378E-1],
                 [1.46820E+2,7.63555E-1],
                 [1.83342E+2,7.47070E-1],
                 [2.40104E+2,7.27294E-1],
                 [3.04907E+2,7.02429E-1],
                 [3.99133E+2,6.86013E-1],
                 [5.14002E+2,6.74616E-1],
                 [6.72413E+2,6.63242E-1],
                 [9.06564E+2,6.61998E-1],
                 [1.04406E+3,6.63885E-1]])
xfr=ipts[:,0]
yfr=ipts[:,1]
frfunc=intp.interp1d(np.log10(xfr),yfr,'cubic')
def fracs(lam):
    """fraction of s poplarization at the ALS
    
    Parameters
    ----------
    lam : number
        wavelength in nm
        
    Returns
    -------
    fr : number
        fraction of s polarization at this wavelength
        
    Example
    -------
    >>> fracs(25.1)
    0.9272135287221943
    """
    ev=np.log10(1239.8/lam)
    return (frfunc(ev)+1)/2

import pandas
from scipy.interpolate import interp1d

class Index:
    """ Index of refraction from volta
    
    Constructor Parameters
    ----------------------
    material : string
        material from base of .nk file on volta
    
    Method
    ------
    at(wavelength) : returns interpolated index at the given
        wavelength in nm. Wavelength can by an np.array
    
    Example
    -------
    alndx=Index('Al')
    lam=np.linspace(10,400,200)
    ndx=alndx.at(lam)
    """
    def __init__(self, material):
        df=pandas.read_table("http://volta.byu.edu/nk/"+material+".nk",
                           comment=';',delim_whitespace=True)
        val=df.values
        lam=val[:,0]/10
        ndx=val[:,1]+val[:,2]*1j
        self.at=interp1d(lam,ndx,'cubic')

from csv import reader

class Log:
    """ Information from the log file for an ALS run
    
    Constructor Parameters
    ----------------------
    path : string
        Directory where the log file is located
    prefix : string
        prefix for the log file name
        
    Attribute
    -------
    Each attribute is a list indexed by run number.
    filename: filename for the data file
    fullname: filename with full path for the data file
    comment: comment for run
    date: date of run
    time: time of run
    gain: log_10 of the gain of the run
    wavelength: wavelength (assumed nm)
    """
    def __init__(self, path:str = 'X:/ALSData/2018/Feb2018',
                 prefix:str = 'Feb2018'):
        self.path = path
        self.filename = ['']
        self.fullname = ['']
        self.comment = ['']
        self.date = ['']
        self.time = ['']
        self.gain = [0]
        self.wavelength = [0]
        mypath = path+'/'+prefix+'.log'
        with open(mypath, newline='') as f:
            f.readline(); # skip first header line
            f.readline(); # skip second header line
            rdr = reader(f, delimiter='\t')
            for row in rdr:
                self.filename.append(row[2]);
                self.fullname.append(path+'/'+row[2])
                self.comment.append(row[3])
                self.date.append(row[4])
                self.time.append(row[5])
                self.gain.append(float(row[8]))
                self.wavelength.append(float(row[26]))
        self.gain = np.array(self.gain)
        self.wavelength = np.array(self.wavelength)

class Run:
    """ The parsed data form a raw data file
    """
    def __init__(self, log: Log, run: int):
        self.var = []
        self.diode = []
        self.m3 = []
        self.beam = []
        with open(log.fullname[run]) as f:
            f.readline();
            rdr = reader(f, delimiter='\t')
            for row in rdr:
                self.var.append(float(row[0]))
                self.diode.append(float(row[1]))
                self.m3.append(float(row[2]))
                self.beam.append(float(row[3]))
        self.var = np.array(self.var)
        self.diode = np.array(self.diode)
        self.m3 = np.array(self.m3)
        self.beam = np.array(self.beam)

class Runs:
    """ A cached collection of Run objects
    
    Constructor Parameters
    ----------------------
    log : Log
        The Log object for this run. By default, it is assumed
        that the data files are stored in the same directory
        as the log file.
    prefix : the prefix for the data files in the runs. A full
        data file name will be <prefix><run number>.dat
        
    Method
    ------
    [] : overloaded index operator returning a Run object. The
        first time a run is referenced, the information is read
        from the raw data. Subsequent accesses used cached data
        from a dictionary.
        
    Example
    -------
    log = Log('X:/ALSData/2018/Feb2018','Feb2018')
    runs = Runs(log)
    run21 = runs[21]
    """
    def __init__(self, log:Log):
        self.log = log
        self.rn = {}
    def __getitem__(self,index):
        sndx = str(index)
        if sndx not in self.rn:
            self.rn[sndx] = Run(self.log, index)
        return self.rn[sndx]