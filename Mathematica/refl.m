(* ::Package:: *)

BeginPackage["refl`"]

refl::usage = "refl is a collection of routines related to finding optical and physical \
properties of thin film optics.  It includes the routines parrattR, parrattRough, matR, \
imdFile, imdMenu, imdIndex, imdIndexAt, imdIndexFunction."

matR::usage = "matR[nd, theta, lambda] returns a two element list with the reflectance\
and transmittance of a stack of thin films on a substrate (which could be a vacuum). \
[At this point, the transmittance is not yet implemented.] The list nd is a list of index \
of refraction/thickness pairs with the substrate as the first element and the vacuum usually\
as the last element. The thickness of the substrate and vacuum layers should be set to zero. \
The thickness should be in nanometers (nm). \
The argument theta is the angle of incidence and reflection measured from grazing on the \
vacuum side in degrees. The lambda argument is the wavelength in nm. \
This calculation is slower than using the Parratt formula if only reflection is desired."

parrattR::usage = "parrattR[nd, theta, lambda] returns the reflectance of a stack of thin \
films on a substrate (which could be a vacuum). The list nd is a list of index of \
refraction/thickness pairs with the substrate as the first element and the vacuum usually \
as the last element.  The thickness of the substrate and vacuum layers should be set to zero. \
The thickness should be in nanometers (nm). \
The argument theta is the angle of incidence and reflection measured from grazing on the \
vacuum side in degrees. The lambda argument is the wavelength in nm. \
This calculation is faster than using the matrix formula in matR if only reflection is desired."

parrattRough::usage = "parrattRough[nd, theta, lambda, sigma] returns the reflectance of a \
stack of thin films on a substrate (which could be a vacuum). \
The list nd is a list of index of refraction/thickness pairs with the substrate as the \
first element and the vacuum usually as the last element. \
The thickness of the substrate and vacuum layers should be set to zero. \
The thickness should be in nanometers (nm). \
The argument theta is the angle of incidence and reflection measured from grazing \
on the vacuum side in degrees. The lambda argument is the wavelength in nm. \
The calculation is corrected for roughness between each layer having a surface rms \
roughness of sigma \
using the Debye-Waller correction for the top layer and the Nevot-Croce correction
between subsequent layers. The sigma list should have one fewer element than nd. \
This calculation is faster than using the matrix formula in matR if only reflection is desired."

imdFile::usage = "imdFile[] returns a list of strings of filenames with index of refraction in them."

imdMenu::usage = "imdMenu[selection] draws a menu button which returns the selected file \
name with index of refraction data as the value of the argument selection on return."

imdIndex::usage = "imdIndex[file] returns a list of the complex index of refraction stored in file."

imdIndexAt::usage = "imdIndexAt[file, wavelength] returns the complex index of refraction from \
file interpolated to the given wavelength measured in nanometers."

imdIndexFunction::usage = "imdIndexFunction[file] returns the interpolation function of the \
index of refraction in file."

Begin["`Private`"]

ddat={{1.05685*^1,0.937888},
{1.40446*^1,0.926536},
{1.75270*^1,0.915094},
{2.25761*^1,0.902016},
{2.95528*^1,0.885600},
{3.80743*^1,0.870842},
{4.67841*^1,0.857696},
{5.74863*^1,0.844550},
{7.06822*^1,0.826362},
{8.69072*^1,0.808173},
{1.15690*^2,0.783378},
{1.46820*^2,0.763555},
{1.83342*^2,0.747070},
{2.40104*^2,0.727294},
{3.04907*^2,0.702429},
{3.99133*^2,0.686013},
{5.14002*^2,0.674616},
{6.72413*^2,0.663242},
{9.06564*^2,0.661998},
{1.04406*^3,0.663885}};
{energy,pol}=Transpose[ddat];
ldat=Transpose[{Log[energy],pol}];

fracS[lambda_]:=Module[{eV, pol, intFunc},
eV=1239.8/lambda; (* Divide hc by lambda to get eV, assuming lambda is in nm *)
intFunc = Interpolation[ldat];
pol=intFunc[Log[eV]];
(pol+1)/2]

matR[nd_, theta_, lambda_]:=Module[{rs,ts, rp, tp, ky,ki, kzi, Ci, A, B,rmats, rm, rmatp,
percentS},
ky=2*Pi*Cos[theta*Pi/180]/lambda;
       ki[ni_]:=2*Pi*ni/lambda;
       kzi[ni_]:=Sqrt[ki[ni]^2-ky^2];
       Ci[ni_, di_]:=Exp[I*kzi[ni]*di/2];
rmats[n1_, n2_, d1_, d2_]:= Module[{kz1, kz2, fs12, gs12, fs21, gs21, C1, C2},
kz1=kzi[n1];
kz2=kzi[n2];
fs12=(kz1-kz2)/(kz1+kz2);
gs12=2kz1/(kz1+kz2);
fs21=-fs12;
gs21=2kz2/(kz1+kz2);
C1=Ci[n1,d1];
C2=Ci[n2 ,d2];
{{gs21*C1*C2-fs21*fs12*C1*C2/gs12,
fs12*C1/(gs12*C2)},
{-fs21*C2/(gs12*C1),
1/(gs12*C1*C2)}}
];
rmatp[n1_, n2_, d1_, d2_]:= Module[{kz1, kz2, fp12, gp12, fp21, gp21, C1, C2},
kz1=kzi[n1];
kz2=kzi[n2];
fp12=(n2^2*kz1-n1^2*kz2)/(n1^2*kz2+n2^2*kz1);
(* Print[StringForm["fp12=``",fp12]]; *)
gp12=2n1*n2*kz1/(n1^2*kz2+n2^2*kz1);
fp21=-fp12;
gp21=2*n1*n2*kz2/(n1^2*kz2+n2^2*kz1);
C1=Ci[n1,d1];
C2=Ci[n2, d2];
{{gp21*C1*C2-fp21*fp12*C1*C2/gp12,
fp12*C1/(gp12*C2)},
{-fp21*C2/(gp12*C1),
1/(gp12*C1*C2)}}
];
A={{1,0},{0,1}};
B=A;
For[i=2, i<=Length[nd], i++,
A=rmats[nd[[i,1]], nd[[i-1,1]],nd[[i,2]],nd[[i-1,2]]].A;
B=rmatp[nd[[i,1]], nd[[i-1,1]],nd[[i,2]],nd[[i-1,2]]].B;
];
ts=1/A[[2,2]];
rs=ts*A[[1,2]];
tp=1/B[[2,2]];
rp=tp*B[[1,2]];
(* Returning the transmission is a bit subtle.  I'll save that for later *)
percentS=fracS[lambda];
percentS*Abs[rs]^2+(1-percentS)*Abs[rp]^2
]

parrattR[nd_, theta_, lambda_]:=Module[{rs, rp, Si, fs, fp, k, Ci, f, Ci12, percentS},
k:=2*Pi/lambda;
Si[ni_]:=Sqrt[ni^2-Cos[theta*Pi/180]^2];
Ci[ni_, di_]:=Exp[2*I*k*Si[ni]*di];
fs[n1_, n2_]:=(Si[n1]-Si[n2])/(Si[n1]+Si[n2]);
fp[n1_, n2_]:=(n2^2*Si[n1]-n1^2*Si[n2])/(n2^2*Si[n1]+n1^2*Si[n2]);
rs=0;
rp=0;
For[i=2, i<=Length[nd],i++,
	f=fs[nd[[i,1]],nd[[i-1,1]]];
	Ci12=Ci[nd[[i,1]],nd[[i,2]]];
	rs=Ci12*(f+rs)/(1+f*rs); 
	f=fp[nd[[i,1]], nd[[i-1,1]]];
	rp=Ci12*(f+rp)/(1+f*rp);
];
percentS=fracS[lambda];
percentS*Abs[rs]^2+(1-percentS)*Abs[rp]^2
]

parrattRough[nd_, theta_, lambda_, sigma_]:=Module[{rs, rp, Si, fs, fp, k, Ci, f, Ci12,
	q1, q2, dbf, percentS, atten},
k:=2*Pi/lambda;
Si[ni_]:=Sqrt[ni^2-Cos[theta*Pi/180]^2];
Ci[ni_, di_]:=Exp[2*I*k*Si[ni]*di];
fs[n1_, n2_]:=(Si[n1]-Si[n2])/(Si[n1]+Si[n2]);
fp[n1_, n2_]:=(n2^2*Si[n1]-n1^2*Si[n2])/(n2^2*Si[n1]+n1^2*Si[n2]);
rs=0;
rp=0;
For[i=2, i<=Length[nd],i++,
	f=fs[nd[[i,1]],nd[[i-1,1]]];
	q1=4 Pi Re[nd[[i-1,1]]]/lambda Sin[theta*Pi/180];
	q2=4 Pi Re[nd[[i,1]]]/lambda Sin[theta*Pi/180];
	atten=Exp[-q1 q2 sigma[[i-1]]^2/2];
	f=f*atten;
	Ci12=Ci[nd[[i,1]],nd[[i,2]]];
	rs=Ci12*(f+rs)/(1+f*rs); 
	f=fp[nd[[i,1]], nd[[i-1,1]]];
	f=f*atten;
	rp=Ci12*(f+rp)/(1+f*rp);
];
percentS=fracS[lambda];
percentS*Abs[rs]^2+(1-percentS)*Abs[rp]^2
]

imdFile[]:=Module[{nkLines,nkFiles},
nkFiles=Import["http://volta.byu.edu/nk","HTML"];
nkLines=StringSplit[nkFiles,"\n"];
First[StringSplit[#]]&/@ nkLines[[6;;-3]]
]

(* Read the index of refraction data from a file *)
dataList[file_String]:=Module[{nkDir="http://volta.byu.edu/nk",nkData},
nkData=Import[StringJoin[{nkDir,"/",file}],"Table"];
(* just keep the header information *)
nkData=Select[nkData,MatchQ[First[#],_Real]&];
(* convert to nm and complex numbers *)
Map[{#[[1]]/10,#[[2]]+I #[[3]]}&,nkData]
]

imdMenu[selection_]:=Module[{names},
names=imdFile[];
PopupMenu[Dynamic[selection],names]
]

imdIndex[file_String]:=Module[{data},
data=dataList[file];
Fold[If[MemberQ[Map[First,#1], First[#2]],#1,Append[#1,#2]]&,{},data]
]

imdIndexFunction[file_String]:=Module[{goodList},
goodList=imdIndex[file];
Interpolation[goodList]
]

imdIndexAt[file_String, wavelength_Real]:=Module[{ifunc},
ifunc=imdIndexFunction[file];
ifunc[wavelength]
]

End[]

EndPackage[ ]
