(* ::Package:: *)

BeginPackage["refl`"]

refl::usage = "refl is a collection of routines related to finding optical and physical properties of thin film optics.  It includes the routines parrattR, parrattRough, matR, imdFile, imdMenu, imdIndex, imdIndexAt, imdIndexFunction."

matR::usage = "matR[nd, theta, lambda, fractionS] returns a two element list with the reflectance and transmittance of a stack of thin films on a substrate (which could be a vacuum).  [At this point, the transmittance is not yet implemented.]  The list nd is a list of index of refraction/thickness pairs with the substrate as the first element and the vacuum usually as the last element.  The thickness of the substrate and vacuum layers should be set to zero.  The argument theta is the angle of incidence and reflection measured from grazing on the vacuum side in degrees.  The lambda argument is the wavelength wavelength using the same units as the layer thickness.  The fractionS argument is the intensity of the s-polarized component of the incident light.  This calculation is slower than using the Parratt formula if only reflection is desired."

parrattR::usage = "parrattR[nd, theta, lambda, fractionS] returns the reflectance of a stack of thin films on a substrate (which could be a vacuum). The list nd is a list of index of refraction/thickness pairs with the substrate as the first element and the vacuum usually as the last element.  The thickness of the substrate and vacuum layers should be set to zero.  The argument theta is the angle of incidence and reflection measured from grazing on the vacuum side in degrees.  The lambda argument is the wavelength wavelength using the same units as the layer thickness.  The fractionS argument is the intensity of the s-polarized component of the incident light.  This calculation is faster than using the matrix formula in matR if only reflection is desired."

parrattRough::usage = "parrattRough[nd, theta, lambda, fractionS, sigma] returns the reflectance of a stack of thin films on a substrate (which could be a vacuum). The list nd is a list of index of refraction/thickness pairs with the substrate as the first element and the vacuum usually as the last element.  The thickness of the substrate and vacuum layers should be set to zero.  The argument theta is the angle of incidence and reflection measured from grazing on the vacuum side in degrees.  The lambda argument is the wavelength wavelength using the same units as the layer thickness.  The fractionS argument is the intensity of the s-polarized component of the incident light.  The calculation is corrected for roughness of the top interface between the vacuum and top layer having a surface rms roughness of sigma using the Debye-Waller correction. This calculation is faster than using the matrix formula in matR if only reflection is desired."

imdFile::usage = "imdFile[] returns a list of strings of filenames with index of refraction in them."

imdMenu::usage = "imdMenu[selection] draws a menu button which returns the selected file name with index of refraction data as the value of the argument selection on return."

imdIndex::usage = "imdIndex[file] returns a list of the complex index of refraction stored in file."

imdIndexAt::usage = "imdIndexAt[file, wavelength] returns the complex index of refraction from file interpolated to the given wavelength measured in nanometers."

imdIndexFunction::usage = "imdIndexFunction[file] returns the interpolation function of the index of refraction in file."

Begin["`Private`"]

matR[nd_, theta_, lambda_, percentS_]:=Module[{rs,ts, rp, tp, ky,ki, kzi, Ci, A, B,rmats, rm, rmatp},
ky=2*Pi*Cos[theta*Pi/180]/lambda;
       ki[ni_]:=2*Pi*ni/lambda;
       kzi[ni_]:=Sqrt[ki[ni]^2-ky^2];
       Ci[ni_, di_]:=Exp[I*kzi[ni]*di/2];
rmats[n1_, n2_, d1_, d2_]:= Module[{kz1, kz2, fs12, gs12, fs21, gs21, C1, C2},
kz1=kzi[n1];
kz2=kzi[n2];
fs12=(kz1-kz2)/(kz1+kz2);
(* Print[StringForm["fs12=``",fs12]]; *)
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
(* Print[A]; *)
For[i=2, i<=Length[nd], i++,
A=rmats[nd[[i,1]], nd[[i-1,1]],nd[[i,2]],nd[[i-1,2]]].A;
B=rmatp[nd[[i,1]], nd[[i-1,1]],nd[[i,2]],nd[[i-1,2]]].B;
];
ts=1/A[[2,2]];
rs=ts*A[[1,2]];
tp=1/B[[2,2]];
rp=tp*B[[1,2]];
(* Returning the transmission is a bit subtle.  I'll save that for later *)
percentS*Abs[rs]^2+(1-percentS)*Abs[rp]^2
]

parrattR[nd_, theta_, lambda_, percentS_]:=Module[{rs, rp, Si, fs, fp, k, Ci, f, Ci12},
k:=2*Pi/lambda;
Si[ni_]:=Sqrt[ni^2-Cos[theta*Pi/180]^2];
Ci[ni_, di_]:=Exp[2*I*k*Si[ni]*di];fs[n1_, n2_]:=(Si[n1]-Si[n2])/(Si[n1]+Si[n2]);
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
(* Print[StringForm["rs=`` rp=``",rs, rp]]; *)
percentS*Abs[rs]^2+(1-percentS)*Abs[rp]^2
]

parrattRough[nd_, theta_, lambda_, percentS_, sigma_]:=Module[{rs, rp, Si, fs, fp, k, Ci, f, Ci12, q, dbf},
k:=2*Pi/lambda;
Si[ni_]:=Sqrt[ni^2-Cos[theta*Pi/180]^2];
Ci[ni_, di_]:=Exp[2*I*k*Si[ni]*di];fs[n1_, n2_]:=(Si[n1]-Si[n2])/(Si[n1]+Si[n2]);
fp[n1_, n2_]:=(n2^2*Si[n1]-n1^2*Si[n2])/(n2^2*Si[n1]+n1^2*Si[n2]);
rs=0;
rp=0;
For[i=2, i<=Length[nd],i++,
f=fs[nd[[i,1]],nd[[i-1,1]]];
If[i==Length[nd],
	q=4 Pi Re[nd[[i-1,1]]]/lambda Sin[theta*Pi/180];
	f=f*Exp[-q^2 sigma^2/2];
];
Ci12=Ci[nd[[i,1]],nd[[i,2]]];
rs=Ci12*(f+rs)/(1+f*rs); 
f=fp[nd[[i,1]], nd[[i-1,1]]];
If[i==Length[nd],
	(* Print["f before ",f]; *)
	f=f*Exp[-q^2 sigma^2/2];
	(* Print["f after ", f]; *)
]; (* sigma^2/2?? *)
rp=Ci12*(f+rp)/(1+f*rp);
];
(* Print[StringForm["rs=`` rp=``",rs, rp]]; *)
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
