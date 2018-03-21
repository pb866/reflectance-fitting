(* ::Package:: *)

BeginPackage["plotAndFit`"]
(* PACKAGE AND FUNCTION USAGE *)


plotAndFit::usage = "Set of routines to plot and fit reflectance as a function \
of angle; contains functions ReflectanceData, FitThicknesses, FitIndex \
(pending), and PlotData

'data' is well-used by this package; it is simply a list of lists, of the form \
{{angle1, reflectance1}, {angle2, reflectance2}, {...}, ...}

This package relies on refl.m; refl` must be imported for it to work properly!"


ReflectanceData::usage = "Take a set of filenames, and optionally a set of \
corresponding angles, and integrate the data therein (subtracting the average \
background)

ReflectanceData[fileNames, angles=None]
    fileNames: a list of the names of the files containing reflectance data
    angles: a list of the angles, from grazing (in degrees), at which \
measurements were taken; if no argument is provided, angles are inferred to be \
the average angular value in the respective files.
    returns data of the form {{angle, reflectance}, {...}, ...}

This function assumes that data are stored as csv/tsv/etc., and are of the form:
    angle 1    refl 1
    angle 2    refl 2
    angle 3    refl 3
    ...
It is thus of limited utility if the file format is not such"


FitThicknesses::usage = "Fit data given substrate material, layer materials \
and estimated thicknesses, wavelength, fraction of S-polarized light, and and \
optionally a weight function

FitThicknesses[data, substrate=\"SiO2\", materialsAndThicknessGuesses, \
weightFunction=\!\(\*FractionBox[\(1\), \(#\)]\)&, \[Lambda]=30.4, fracS=0.5]
    data: reflectance of angle (degrees); stored as {{angle,refl}, {...}, ...}
    substrate: substrate material; string, e.g. \"SiO2\"; silicon by default
    materialsAndThicknessGuesses: of the form {\"Si3N4\", 300, SiNThickness}, \
{\"Al\", 50, AlThickness}} (from substrate to surface; vacuum automatically \
appended); thickness in nm, and SiNThicnkess, etc. are variables for the fit
    weightFunction: weighting function; how much weight should be given to a \
fit parameter as a function of reflectance?
    \[Lambda]: wavelength (in nm) of incident light; default 30.4
    fracS: how much of light is s-polarized? Default 0.5
    returns a NonlinearModelFit of the data, of the form of the Parratt \
recursion function (which must be defined by importing refl`"


PlotData::usage = "ListLogPlot data and, optionally, a fit

ListLogPlot[data, fit=None]
    displays a log plot of the data points and, if specified, the fit"


(* THE FUNCTIONS THEMSELVES *)
Begin["`Private`"]


IntegrateData[data_, incidentAngle_:None, wingSize_:1] := Module[
    {leftWing,rightWing,leftAvg,rightAvg,middle,middleT,background,vals,angle},
    (*
    data: data from Imported file to be integrated
    wingSize: how much of edges do we want to average over to get out noise? (in degrees; default .5)
    returns a pair--{angle,sum}
    *)

    (*Find the wings and use them to find the average background*)
    leftWing=Select[data,First[data][[1]]-#[[1]]+wingSize>0&];
    rightWing=Select[data,#[[1]]-Last[data][[1]]+wingSize>0&];
    leftAvg=Mean[Transpose[leftWing][[2]]];
    rightAvg=Mean[Transpose[rightWing][[2]]];
    background=(leftAvg+rightAvg)/2//N;
    
    (*Use the middle portion, subtracting the background, to "integrate" the data*)
    middle=DeleteCases[data,Alternatives@@Join[leftWing,rightWing]];
    middleT=Transpose[middle];
    angle=-Mean[middleT[[1]]]/2;
    
    (*Subtract background from middle*)
    vals=#-background&/@(middleT[[2]]);
    
    (*Return the angle and total of the middle portion*)
    {If[incidentAngle==None,angle,incidentAngle],Total[vals]}
]


ReflectanceData[fileNames_,angles_:None]:=Module[{data,dataT,refls},
    (*
    fileNames: a list of file names in which data are found
    angles: angles of samples at measurement
    returns data normalized to 1
    *)
    
    (*Import and integrate data from a given file*)
    data=Import[#,"Data"]&/@fileNames;
    data=If[angles==None,
        IntegrateData[#]&/@data,
        (*ELSE*)
        IntegrateData[#[[1]],#[[2]]]&/@{data,angles}//Transpose
    ];
    dataT=data//Transpose;
    
    (*Extract reflection column*)
    refls=#/Max[dataT[[2]]]&/@dataT[[2]];
    
    (*`data` is now of the format {{angle,reflectance},{...}...} where angle is angle FROM GRAZING*)
    data={dataT[[1]],refls}//Transpose;
    
    (*Discard any reflections less than zero (since they're likely garbage) and return*)
    Select[data,#[[2]]>0&]
]


FitThicknesses[data_,substrate_,materialsAndThicknessGuesses_:{},weightFunction_: 1/#&,\[Lambda]_:30.4,fracS_:0.5]:=Module[{materialFiles,thicknesses,thicknessVars,indices,layers,weights,constraints,thicknessesAndGuesses},
    (*
    data: reflectance of angle (in degrees); stored as {{angle,refl},{...},...}
    substrate: the material of the substrate; string, e.g. "SiO2"
    materialsAndThicknessGuesses: of the form {"Si3N4",0,SiNThickness},{"Al",50,AlThickness}} (from substrate to surface; vacuum automatically appended, substrate thickness guess ignored); thickness in nm, and SiNThicnkess etc. are variables for the fit
    weightFunction: weighting function; how much weight should be given to a fit parameter as a function of reflectance?
    \[Lambda]: wavelength (in nm) of incident light for this sample; default 30.4
    fracS: how much of light is s-polarized? Default 0.5
    *)
    
    (*Get material filenames TODO: initialize, setting notebook directory to this directory and importing refl*)
    materialFiles={substrate<>".nk"}~Join~(#[[1]]<>".nk"&/@materialsAndThicknessGuesses);
    
    (*Get thicknesses and indices*)
    thicknesses={0}~Join~(#[[2]]&/@materialsAndThicknessGuesses);
    
    (*Drop first thickness guess for our attempt to fit--we don't need it*)
    thicknessVars=#[[3]]&/@materialsAndThicknessGuesses;
    indices=imdIndexAt[#,\[Lambda]]&/@materialFiles;
    
    (*Combine thicknesses and indices into layers; format is {{index,thicknessVariable},{...},...}; also appends a vacuum layer*)
    layers=({indices,{0}~Join~thicknessVars}//Transpose)~Join~{{1,0}};
    
    (*Compute weights*)
    weights=weightFunction&/@(data//Transpose)[[2]];
    
    (*Constraints (just making sure thicknesses are greater than zero)*)
    constraints=thicknessVars[[1]]>0;
    For[i=2,i<=Length[thicknessVars],i++,
        constraints=constraints&&thicknessVars[[i]]>0
    ];
    
    (*Thickness variables, and our guesses for their thicknesses*)
    thicknessesAndGuesses={thicknessVars,Drop[thicknesses,1]}//Transpose;
    
    (*Return a fit*)
    NonlinearModelFit[data,{parrattR[layers,\[Theta],\[Lambda],fracS],constraints},thicknessesAndGuesses,\[Theta],Weights->weights]
]


PlotData[data_,fit_:None]:=Module[{MyListPlot,layers,i,thickness,\[Theta]},
    (*
    data: reflectance of angle (in degrees); stored as {{angle,refl},{...},...}
    wavelength: wavelength (in nm) of incident light
    fracS: amount of incident light which is s-polarized
    layerMaterials: list of layer materials (from substrate to surface; vacuum is automatically appended) constituting the thin film; defaults to Null, in which case no fit is attempted
    *)
    MyListPlot:=ListLogPlot[data,PlotLabel->"Reflectance vs Angle from Grazing",AxesLabel->{"Angle (deg)","Reflectance (arb)"},ImageSize->Large];
    Show[MyListPlot,LogPlot[fit[\[Theta]]//Normal,{\[Theta],0,90}]]
]


End[] (*Private*)
EndPackage[]
