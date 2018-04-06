(* ::Package:: *)

BeginPackage["rdata`"]
rdata::usage="rdata is a collection of routines for combining ALS run data to generate reflectance data. For now, it has the single routine reflData."
reflData::usage="reflData[iRuns,i0run,lambda,darkRuns,dropPoints] Combine run data from ALS to make a reflectance file. iRuns is a list of run numbers for I, i0run is the run number with the I0 data, lambda is the wavelength in nm, dark runs are the run numbers for the dark current (one run number for each gain), dropPoints are the number of points to drop at the beginning of the spectrum. The function will return a sorted list of pairs of points with reflectance as a function of angle."
Begin["`Private`"]
path="..\\Y2O3\\1504y2o3 (day 1)\\";
prefix="y2o3001";
postfix=".dat";
fileName[run_]:=path<>prefix<>IntegerString[run,10,3]<>postfix;
logFileName=path<>"1504y2o3.log";

runGain[run_]:=Module[{logData, gain},
logData=Import[logFileName];
gain=Drop[logData,2][[run+1, 9]];
10^gain
]

readRun[run_,dark_]:=Module[{theta, i, allData, current, gain,r},
allData=Drop[Import[fileName[run],"Table"],1]//Transpose;
allData = Transpose[Select[Transpose[allData], #[[2]] < 10 &]];
theta=allData[[1]];
i=allData[[2]];
current=allData[[4]];
gain=runGain[run];
r=(i-dark)/gain/current;
Transpose[{theta,r}]
]

computeDark[darkRuns_, gain_]:=Module[{gainRun,darkData},
gainRun=Select[darkRuns,runGain[#]==gain&];
darkData=Drop[Import[fileName[gainRun],"Table"],1]//Transpose;
Mean[darkData[[2]]]
]

reflData[iRuns_,i0run_,lambda_,darkRuns_,dropPoints_]:=Module[{iData,i0Func, i0,rData, tData}, 
iData=Catenate[readRun[#,computeDark[darkRuns,runGain[#]]]&/@iRuns];
i0Func=Interpolation[readRun[i0run,computeDark[darkRuns, runGain[i0run]]]];
i0=i0Func[lambda];
tData=Transpose[iData];
rData=tData[[2]]/i0;
iData=Drop[Sort[Transpose[{tData[[1]],rData}]],dropPoints] 
]

End[]
EndPackage[ ]
