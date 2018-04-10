(* ::Package:: *)

(* ::Title:: *)
(*Sample 7 Analysis*)


(* ::Chapter:: *)
(*Gabriel Richardson*)


(* ::Text:: *)
(**)


(* ::Section:: *)
(*Setup*)


(* ::Text:: *)
(*First we need to set the default directory. I'll also show how to develop a small module to read in a data file and return the columns of data stored in it.*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*<<refl`*)


(* ::Input:: *)
(*"/Users/GabrielRichardson/Desktop/Research/ALS Analysis"*)


(* ::Text:: *)
(*The following module reads in a data file specified by run number and returns the four columns of data in the file. The prefix is the location of the files and the start of the file name. The postfix is the extension. These are concatenated together using the <> operators together with the run number which is passed to the module. It needs to be converted to a string before being concatenated. The Import command reads in the file as a Table of rows. The Drop command drops the first row of the table which is a comment with the date and time. The Transpose commands converts the data from a list of rows to a list of columns which is what the module returns.*)


(* ::Input:: *)
(*readRun[run_]:=Module[{prefix="alf3000",postfix=".txt",data},*)
(*data=Drop[Import[prefix<>ToString[run]<>postfix,"Table"],1];*)
(*Transpose[data]*)
(*]*)


(* ::Text:: *)
(*Here is an example of its use:*)


(* ::Input:: *)
(*{c1,c2,c3,c4}=readRun[366];*)


(* ::Text:: *)
(*Run 359 is the I0 scan, so c1 is the wavelength and c2 the I0 signal. I can use LIstPlot to plot this. I'll need to take the two columns I need, combine them, and then Transpose them back to a list of x,y pairs of data.*)


(* ::Input:: *)
(*ListPlot[Transpose[{c1,c2}]]*)


(* ::Input:: *)
(*ListPlot[Transpose[{c1,c2}],Joined->True, PlotLabel->"I0",AxesLabel->{"\[Lambda], nm","current"}]*)


(* ::Section:: *)
(*Dark Current*)


(* ::Input:: *)
(**)
(*(*You can change the runs to correspond with the dark current runs you are using, mine are 367-369*)*)


(* ::Input:: *)
(*{c1,c2,c3,c4}=readRun[368];*)
(*d9=Mean[c2]*)
(*{c1,c2,c3,c4}=readRun[367];*)
(*d8=Mean[c2]*)
(*{c1,c2,c3,c4}=readRun[369];*)
(*d10=Mean[c2]*)


(* ::Section:: *)
(*I0*)


(* ::Input:: *)
(*(*My io run was 366, but you could change it to yours*)*)
(**)
(*{\[Lambda],inaught,c3,c4}=readRun[366];*)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[17.1]*)


(* ::Input:: *)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[17.1]*)


(* ::Input:: *)
(*7.342864703752149`*)


(* ::Section:: *)
(*Index of Refraction*)


(* ::Input:: *)
(*(*For any other element for which you may need the index of refraction, you can follow the same model as is for the others bellow*)*)
(**)
(*sindx=0.8547067+0.01427I*)
(*sio2ndx=0.86515739+0.19722387I*)
(**)
(*si=Import["Si.nk","Table"];*)
(*si=Drop[si,8];*)
(*si[[1;;4]]*)
(*{wl,n,k}=Transpose[si];*)
(*si=Transpose[{wl,n+k*I}];*)
(*si[[1;;4]]*)
(*siInterp=Interpolation[si];*)
(**)
(*siO2=Import["SiO2.nk","Table"];*)
(*siO2=Drop[siO2,8];*)
(*siO2[[1;;4]]*)
(*{wl,n,k}=Transpose[siO2];*)
(*siO2=Transpose[{wl,n+k*I}];*)
(*siO2Interp=Interpolation[siO2];*)
(**)
(*alF3=Import["AlF3.nk","Table"];*)
(*alF3=Drop[alF3,8];*)
(*{wl,n,k}=Transpose[alF3];*)
(*alF3=Transpose[{wl,n+k*I}];*)
(*alF3Interp=Interpolation[alF3];*)
(**)
(*si3n4=Import["Si3N4.nk","Table"];*)
(*si3n4=Drop[si3n4,8];*)
(*{wl,n,k}=Transpose[si3n4];*)
(*si3n4=Transpose[{wl,n+k*I}];*)
(*si3n4Interp=Interpolation[si3n4,InterpolationOrder->1];*)
(**)
(*al=Import["Al.nk","Table"];*)
(*al=Drop[al,8];*)
(*{wl,n,k}=Transpose[al];*)
(*al=Transpose[{wl,n+k*I}];*)
(*alInterp=Interpolation[al,InterpolationOrder->1];*)


(* ::Section:: *)
(*Models*)


(* ::Input:: *)
(*(*This will compute the reflectance*)*)
(*compute[run1_,run2_,lda_,gain2_]:=Module[{},*)
(*{\[Theta],iDatan1,c3,c4}=readRun[run1];*)
(**)
(*{\[Theta]2,iDatan2,c3,c4}=readRun[run2]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iDatan1~Join~(iDatan2*10^-(gain2-8)); *)
(**)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[lda];*)
(*refla=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[lda];*)
(*Export["S7W"<>ToString[lda]<>".txt",Transpose[{\[Theta],refla}],"Table"];*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refla}]]*)
(**)
(*]*)
(**)
(*(*This will perform the fitting, but it is created according to this sample. This is where most of the modification will come in*)*)
(*fitting[ang_,th_,ndx_,beta_,theta_]:=Module[{},*)
(*nd={{siInterp[ang],0.0},{siO2Interp[ang],16.0},{si3n4Interp[ang],1090.0},{alInterp[ang],130.0},{alF3Interp[ang],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[ang],0.0},{alInterp[ang],130.0},{alF3Interp[ang],85.0},{1.0,0.0}};*)
(*lam=ang;*)
(*fractionS=0.9;*)
(*r=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*rb=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t, rd}=Transpose[r];*)
(*{t, rdb}=Transpose[rb];*)
(*rdif=Transpose[{t,(rd-rdb)}];*)
(*p=ListLogPlot[{r,rb}, Joined->True];*)
(*g=ListPlot[rdif];*)
(*alF3Interp[ang];*)
(*nd={{si3n4Interp[ang],0.0},{alInterp[ang],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d=Drop[Import["S7W"<>ToString[ang*10^-1]<>".txt","Table"],5];*)
(*{angl,wt}=Transpose[d];*)
(*weights=1/wt^2;*)
(*fit=NonlinearModelFit[d,parrattR[nd,theta,lam,fractionS],{{th,130.0}, {ndx,0.936}, {beta,0.102}}, {theta},Weights->weights];*)
(*y=fit["ParameterTable"];*)
(*z=Show[LogPlot[fit[theta],{theta,0,80}],ListLogPlot[d]];*)
(*CellPrint@ExpressionCell[#,"Output"]&/@{p,g,y,z};*)
(*]*)
(**)


(* ::Section:: *)
(*17.1 nm*)


(* ::Input:: *)
(*{\[Theta],iData,c3,c4}=readRun[370];*)


(* ::Input:: *)
(*refl=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)


(* ::Input:: *)
(*ListLogPlot[Transpose[{\[Theta],refl}]]*)


(* ::Input:: *)
(*Export["S7W17m.txt",Transpose[{\[Theta],refl}],"Table"]*)


(* ::Section:: *)
(*18 nm*)


(* ::Input:: *)
(*{\[Theta]1,iData1,c3,c4}=readRun[372];*)
(**)


(* ::Input:: *)
(*{\[Lambda],inaught,c3,c4}=readRun[366];*)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[18]*)
(*refl18=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[18]*)


(* ::Input:: *)
(*{\[Theta]2,iData2,c3,c4}=readRun[373]; (*This is the same wavelength, just at gain 10*)*)
(*{\[Lambda],inaught,c3,c4}=readRun[366];*)
(**)


(* ::Input:: *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iData1~Join~(iData2*10^-2); (*This makes it so that we can plot both gains*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*ListLogPlot[Transpose[{\[Theta],refl18}]]*)


(* ::Input:: *)
(*Export["S7W18M.txt",Transpose[{\[Theta],refl18}],"Table"]*)


(* ::Input:: *)
(*siInterp[180]*)
(*siO2Interp[180]*)
(*alF3Interp[180]*)
(*si3n4Interp[180]*)
(*alInterp[180];*)
(*(*Can adjust the 16*)*)
(*nd={{siInterp[180],0.0},{siO2Interp[180],18.0},{si3n4Interp[180],1090.0},{alInterp[180],130.0},{alF3Interp[180],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[180],0.0},{alInterp[180],130.0},{alF3Interp[180],85.0},{1.0,0.0}};*)
(*lam=180;*)
(*fractionS=0.95;*)
(*r18=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r18b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t18, r18d}=Transpose[r18];*)
(*{t18, r18db}=Transpose[r18b];*)
(*r18dif=Transpose[{t18,(r18d-r18db)}];*)
(*ListLogPlot[{r18,r18b}, Joined->True]*)


(* ::Input:: *)
(*ListPlot[r18dif]*)


(* ::Input:: *)
(*alF3Interp[180];*)
(*(*Include SiO2 and Si in the nd fit. Si for substrate with SiO2 layer below, just like above.*)*)
(*nd={{si3n4Interp[180],1100.0},{alInterp[180],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d18=Drop[Import["S7W18M.txt","Table"],7];*)
(*{ang,wt}=Transpose[d18];*)
(*weights=1/wt^2;*)
(*(*Need to calculate fractionS with new code. Between .85 and .95*)*)
(*(*fit18=NonlinearModelFit[d18,parrattR[nd,theta,lam,fractionS],{{th,112.0}, {ndx,0.951}, {beta,0.0205}}, {theta}];*)
(*fit18["ParameterTable"]*)*)
(*fit18=NonlinearModelFit[d18,parrattR[nd,theta,lam,fractionS],{{th,112.0}, {ndx,0.951}, {beta,0.0205}}, {theta},Weights->weights];*)
(*fit18["ParameterTable"]*)
(*Show[LogPlot[fit18[theta],{theta,0,80}],ListLogPlot[d18]]*)


(* ::Input:: *)
(*fit18["BestFitParameters"]*)


(* ::Input:: *)
(*{th->104.53615055839997`,ndx->0.9635957796340774`,beta->0.036499299035447415`}*)
(*th18=fit18["BestFitParameters"][[1]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*20 nm*)


(* ::Input:: *)
(*{\[Theta],iData20a,c3,c4}=readRun[374];*)
(*{\[Lambda],inaught,c3,c4}=readRun[366];*)
(**)
(*{\[Theta]2,iData20b,c3,c4}=readRun[375]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iData20a~Join~(iData20b*10^-2); *)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[20]*)
(*refl20=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[20]*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refl20}]]*)


(* ::Input:: *)
(*Export["S7W20M.txt",Transpose[{\[Theta],refl20}],"Table"]*)


(* ::Input:: *)
(*?parrattRough*)


(* ::Input:: *)
(*Clear[ndx]*)
(*siInterp[200];*)
(*siO2Interp[200];*)
(*alF3Interp[200];*)
(*si3n4Interp[200];*)
(*alInterp[200];*)
(*nd={{siInterp[200],0.0},{siO2Interp[200],16.0},{si3n4Interp[200],1090.0},{alInterp[200],130.0},{alF3Interp[200],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[200],0.0},{alInterp[200],130.0},{alF3Interp[200],85.0},{1.0,0.0}};*)
(*lam=200;*)
(*fractionS=0.9;*)
(*r20=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r20b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t20, r20d}=Transpose[r20];*)
(*{t20, r20db}=Transpose[r20b];*)
(*r20dif=Transpose[{t20,(r20d-r20db)}];*)
(*ListLogPlot[{r20,r20b}, Joined->True]*)
(*ListPlot[r20dif]*)
(*alF3Interp[200];*)
(*nd={{siInterp[200],0.0},{siO2Interp[200],16.0},{si3n4Interp[200],1210},{alInterp[200],th},(*{al2O3Interp} Try putting in a layer of oxide to see if it helps, start with just a couple of nm thickness*){ndx+beta*I,85.0},{1.0,0.0}};*)
(*d201=Take[Import["S7W20M.txt","Table"],{8,80}];*)
(*d202=Take[Import["S7W20M.txt","Table"],{86,160}];*)
(*d20=d201~Join~d202;*)
(*{ang,wt}=Transpose[d20];*)
(*weights=1/wt^2;*)
(*(*fit20=NonlinearModelFit[d20,parrattR[nd,theta,lam,fractionS],{{th,130.0}, {ndx,0.936}, {beta,0.102}}, {theta}];*)
(*fit20["ParameterTable"]*)*)
(*fit20=NonlinearModelFit[d20,parrattR[nd,theta,lam,.92],{{th1,1100},{th,75}, {ndx,0.963}, {beta,0.036}}, {theta},Weights->weights];*)
(*fit20["ParameterTable"]*)
(*Show[LogPlot[fit20[theta],{theta,0,80}],ListLogPlot[d20]]*)


(* ::Input:: *)
(*(*Only do where we see wiggles, and see if we can come up with thickness of nitride. Then put the actual substrates for all of them and make the nitride the average of the ones with the wiggles, and we how well we can get the rest fit*)*)
(*(*If this doesnt work, try fitting the index of the aluminum. This could show poor aluminum*)*)


(* ::Input:: *)
(**)
(*th18=th/.fit18["BestFitParameters"][[1]];*)
(*th20=th/.fit20["BestFitParameters"][[1]];*)
(*th22=th/.fit22["BestFitParameters"][[1]];*)
(*th24=th/.fit24["BestFitParameters"][[1]];*)
(*th26=th/.fit26["BestFitParameters"][[1]];*)
(*thicknesses={th18,th20,th22,th24,th26};*)
(*ListPlot[thicknesses]*)
(*Plot[thicknesses,{x,18,26},PlotLegends->"Expressions"]*)


(* ::Input:: *)
(*n18=ndx/.fit18["BestFitParameters"][[2]];*)
(*n20=ndx/.fit20["BestFitParameters"][[2]];*)
(*n22=ndx/.fit22["BestFitParameters"][[2]];*)
(*n24=ndx/.fit24["BestFitParameters"][[2]];*)
(*n26=ndx/.fit26["BestFitParameters"][[2]];*)
(*ns={n18,n20,n22,n24,n26};*)
(*ListPlot[ns]*)


(* ::Input:: *)
(*k18=beta/.fit18["BestFitParameters"][[3]];*)
(*k20=beta/.fit20["BestFitParameters"][[3]];*)
(*k22=beta/.fit22["BestFitParameters"][[3]];*)
(*k24=beta/.fit24["BestFitParameters"][[3]];*)
(*k26=beta/.fit26["BestFitParameters"][[3]];*)
(*ks={k18,k20,k22,k24,k26};*)
(*ListPlot[ks]*)


(* ::Section:: *)
(*22 nm*)


(* ::Input:: *)
(*{\[Theta],iData22a,c3,c4}=readRun[376];*)
(*(*{\[Lambda],inaught,c3,c4}=readRun[366];*)*)
(**)
(*{\[Theta]2,iData22b,c3,c4}=readRun[377]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iData22a~Join~(iData22b*10^-2); *)
(**)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[22]*)
(*refl22=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[22]*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refl22}]]*)


(* ::Input:: *)
(*Export["S7W22.txt",Transpose[{\[Theta],refl20}],"Table"]*)


(* ::Input:: *)
(*siInterp[220];*)
(*siO2Interp[220];*)
(*alF3Interp[220];*)
(*si3n4Interp[220];*)
(*alInterp[220];*)
(*nd={{siInterp[220],0.0},{siO2Interp[220],16.0},{si3n4Interp[220],1090.0},{alInterp[220],130.0},{alF3Interp[220],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[220],0.0},{alInterp[220],130.0},{alF3Interp[220],85.0},{1.0,0.0}};*)
(*lam=220;*)
(*fractionS=0.9;*)
(*r22=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r22b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t22, r22d}=Transpose[r22];*)
(*{t22, r22db}=Transpose[r22b];*)
(*r22dif=Transpose[{t22,(r22d-r22db)}];*)
(*ListLogPlot[{r22,r22b}, Joined->True]*)
(*ListPlot[r22dif]*)
(*alF3Interp[220];*)
(*nd={{siInterp[220],10000.0},{si3n4Interp[220],1100.0},{alInterp[220],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*(*Changig si3n4 from 1090 to 0 we can see that si3n4 is more transparent then we suspected*)*)
(*d22=Drop[Import["S7W22.txt","Table"],5];*)
(*d221=Take[Import["S7W22.txt","Table"],{8,80}];*)
(*d222=Take[Import["S7W22.txt","Table"],{86,160}];*)
(*d22=d221~Join~d222;*)
(*{ang,wt}=Transpose[d22];*)
(*weights=1/wt^2;*)
(*(*fit22=NonlinearModelFit[d22,parrattR[nd,theta,lam,fractionS],{{th,130.0}, {ndx,0.936}, {beta,0.102}}, {theta}];*)
(*fit22["ParameterTable"]*)*)
(*fit22=NonlinearModelFit[d22,parrattR[nd,theta,lam,fractionS],{{th,103.6885}, {ndx,0.9559}, {beta,0.04725}}, {theta},Weights->weights];*)
(*fit22["ParameterTable"]*)
(*Show[LogPlot[fit22[theta],{theta,0,80}],ListLogPlot[d22]]*)


(* ::Input:: *)
(*(*Drop the straight line of data*)*)
(*fit22["BestFitParameters"]*)


(* ::Section:: *)
(*24 nm*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*compute[378,379,24,9]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*siInterp[240];*)
(*siO2Interp[240];*)
(*alF3Interp[240];*)
(*si3n4Interp[240];*)
(*alInterp[240];*)
(*nd={{siInterp[240],0.0},{siO2Interp[240],16.0},{si3n4Interp[240],1090.0},{alInterp[240],130.0},{alF3Interp[240],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[240],0.0},{alInterp[240],130.0},{alF3Interp[240],85.0},{1.0,0.0}};*)
(*lam=240;*)
(*fractionS=0.9;*)
(*r24=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r24b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t24, r24d}=Transpose[r24];*)
(*{t24, r24db}=Transpose[r24b];*)
(*r24dif=Transpose[{t24,(r24d-r24db)}];*)
(*ListLogPlot[{r24,r24b}, Joined->True]*)
(*ListPlot[r24dif]*)
(*alF3Interp[240];*)
(*nd={{si3n4Interp[240],0.0},{alInterp[240],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d24=Drop[Import["S7W24.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d24];*)
(*weights=1/wt1^2;*)
(*fit24=NonlinearModelFit[d24,parrattR[nd,theta,lam,fractionS],{{th,106.51}, {ndx,.9764}, {beta,0.0655}}, {theta},Weights->weights];*)
(*(*Use the guess from the pervious fit for these fits. Ex. beta of .0364, ndx...We can change the module so that it uses the previous fits for the next theta*)*)
(*fit24["ParameterTable"]*)
(*Show[LogPlot[fit24[theta],{theta,0,80}],ListLogPlot[d24]]*)


(* ::Input:: *)
(**)
(*(*The reason it doesnt go all the way down is because there is too much absorption. Get rid of the tailing off at the beginning*)*)


(* ::Section::Closed:: *)
(*26 nm*)


(* ::Input:: *)
(*compute[380,381,26,9]*)


(* ::Input:: *)
(*fit24["BestFitParameters"]*)
(*siInterp[260];*)
(*siO2Interp[260];*)
(*alF3Interp[260];*)
(*si3n4Interp[260];*)
(*alInterp[260];*)
(*nd={{siInterp[260],0.0},{siO2Interp[260],16.0},{si3n4Interp[260],1090.0},{alInterp[260],130.0},{alF3Interp[260],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[260],0.0},{alInterp[260],130.0},{alF3Interp[260],85.0},{1.0,0.0}};*)
(*lam=260;*)
(*fractionS=0.9;*)
(*r26=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r26b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t26, r26d}=Transpose[r26];*)
(*{t26, r26db}=Transpose[r26b];*)
(*r26dif=Transpose[{t26,(r26d-r26db)}];*)
(*ListLogPlot[{r26,r26b}, Joined->True]*)
(*ListPlot[r26dif]*)
(*alF3Interp[260];*)
(*nd={{si3n4Interp[260],0.0},{alInterp[260],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d26=Drop[Import["S7W26.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d26];*)
(*weights=1/wt1^2;*)
(*fit26=NonlinearModelFit[d26,parrattR[nd,theta,lam,fractionS],{{th,101.58}, {ndx,.946465}, {beta,0.0667761}}, {theta},Weights->weights];*)
(*(*I had to remove ,Weights\[Rule]weights because it wasnt working*)*)
(*(*Use the guess from the pervious fit for these fits. Ex. beta of .0364, ndx...We can change the module so that it uses the previous fits for the next theta*)*)
(*fit26["ParameterTable"]*)
(*Show[LogPlot[fit26[theta],{theta,0,80}],ListLogPlot[d26]]*)


(* ::Section::Closed:: *)
(*28*)


(* ::Input:: *)
(*compute[382,383,28,9]*)


(* ::Input:: *)
(*fit26["BestFitParameters"]*)
(*nd={{siInterp[280],0.0},{siO2Interp[280],16.0},{si3n4Interp[280],1090.0},{alInterp[280],130.0},{alF3Interp[280],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[280],0.0},{alInterp[280],130.0},{alF3Interp[280],85.0},{1.0,0.0}};*)
(*lam=280;*)
(*fractionS=0.9;*)
(*r28=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r28b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t28, r28d}=Transpose[r28];*)
(*{t28, r28db}=Transpose[r28b];*)
(*r28dif=Transpose[{t28,(r28d-r28db)}];*)
(*ListLogPlot[{r28,r28b}, Joined->True]*)
(*ListPlot[r28dif]*)
(*nd={{si3n4Interp[280],0.0},{alInterp[280],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d28=Drop[Import["S7W28.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d28];*)
(*weights=1/wt1^2;*)
(*fit28=NonlinearModelFit[d28,parrattR[nd,theta,lam,fractionS],{{th,100.271}, {ndx,.945774}, {beta,0.0590564}}, {theta},Weights->weights];*)
(*fit28["ParameterTable"]*)
(*Show[LogPlot[fit28[theta],{theta,0,80}],ListLogPlot[d28]]*)


(* ::Input:: *)
(**)


(* ::Section::Closed:: *)
(*30*)


(* ::Input:: *)
(*compute[384,385,30,10]*)


(* ::Input:: *)
(*fit28["BestFitParameters"]*)
(*nd={{siInterp[300],0.0},{siO2Interp[300],16.0},{si3n4Interp[300],1090.0},{alInterp[300],130.0},{alF3Interp[300],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[300],0.0},{alInterp[300],130.0},{alF3Interp[300],85.0},{1.0,0.0}};*)
(*lam=300;*)
(*fractionS=0.9;*)
(*r30=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r30b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t30, r30d}=Transpose[r30];*)
(*{t30, r30db}=Transpose[r30b];*)
(*r30dif=Transpose[{t30,(r30d-r30db)}];*)
(*ListLogPlot[{r30,r30b}, Joined->True]*)
(*ListPlot[r30dif]*)
(*nd={{si3n4Interp[300],0.0},{alInterp[300],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d30=Drop[Import["S7W30.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d30];*)
(*weights=1/wt1^2;*)
(*fit30=NonlinearModelFit[d30,parrattR[nd,theta,lam,fractionS],{{th,103.831970983}, {ndx,0.93709435577}, {beta,0.066177}}, {theta},Weights->weights];*)
(*fit30["ParameterTable"]*)
(*Show[LogPlot[fit30[theta],{theta,0,80}],ListLogPlot[d30]]*)


(* ::Section::Closed:: *)
(*New Dark Current and I0 for Higher Wavelengths*)


(* ::Input:: *)
(*{c1,c2,c3,c4}=readRun[388];*)
(*d9=Mean[c2]*)
(*{c1,c2,c3,c4}=readRun[387];*)
(*d8=Mean[c2]*)
(*{c1,c2,c3,c4}=readRun[389];*)
(*d10=Mean[c2]*)
(**)
(*{\[Lambda],inaught,c3,c4}=readRun[386];*)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)


(* ::Section::Closed:: *)
(*26*)


(* ::Input:: *)
(*compute[390,391,26,9]*)


(* ::Input:: *)
(*fit24["BestFitParameters"]*)
(*siInterp[260];*)
(*siO2Interp[260];*)
(*alF3Interp[260];*)
(*si3n4Interp[260];*)
(*alInterp[260];*)
(*nd={{siInterp[260],0.0},{siO2Interp[260],16.0},{si3n4Interp[260],1090.0},{alInterp[260],130.0},{alF3Interp[260],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[260],0.0},{alInterp[260],130.0},{alF3Interp[260],85.0},{1.0,0.0}};*)
(*lam=260;*)
(*fractionS=0.9;*)
(*r26=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r26b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t26, r26d}=Transpose[r26];*)
(*{t26, r26db}=Transpose[r26b];*)
(*r26dif=Transpose[{t26,(r26d-r26db)}];*)
(*ListLogPlot[{r26,r26b}, Joined->True]*)
(*ListPlot[r26dif]*)
(*alF3Interp[260];*)
(*nd={{si3n4Interp[260],0.0},{alInterp[260],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d26=Drop[Import["S7W26.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d26];*)
(*weights=1/wt1^2;*)
(*fit26=NonlinearModelFit[d26,parrattR[nd,theta,lam,fractionS],{{th,101.62}, {ndx,.946407}, {beta,0.0667544}}, {theta},Weights->weights];*)
(*fit26["ParameterTable"]*)
(*Show[LogPlot[fit26[theta],{theta,0,80}],ListLogPlot[d26]]*)


(* ::Section::Closed:: *)
(*28*)


(* ::Input:: *)
(*compute[392,393,28,9]*)


(* ::Input:: *)
(*fit26["BestFitParameters"]*)
(*nd={{siInterp[280],0.0},{siO2Interp[280],16.0},{si3n4Interp[280],1090.0},{alInterp[280],130.0},{alF3Interp[280],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[280],0.0},{alInterp[280],130.0},{alF3Interp[280],85.0},{1.0,0.0}};*)
(*lam=280;*)
(*fractionS=0.9;*)
(*r28=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r28b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t28, r28d}=Transpose[r28];*)
(*{t28, r28db}=Transpose[r28b];*)
(*r28dif=Transpose[{t28,(r28d-r28db)}];*)
(*ListLogPlot[{r28,r28b}, Joined->True]*)
(*ListPlot[r28dif]*)
(*nd={{si3n4Interp[280],0.0},{alInterp[280],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d28=Drop[Import["S7W28.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d28];*)
(*weights=1/wt1^2;*)
(*fit28=NonlinearModelFit[d28,parrattR[nd,theta,lam,fractionS],{{th,100.411}, {ndx,.94544}, {beta,0.0598566}}, {theta},Weights->weights];*)
(*fit28["ParameterTable"]*)
(*Show[LogPlot[fit28[theta],{theta,0,80}],ListLogPlot[d28]]*)
(*(*Something interesting that I noticed as I was inputing the data was that this 28 was not much different from the last, same with 26*)*)


(* ::Section::Closed:: *)
(*30*)


(* ::Input:: *)
(*compute[394,395,30,9]*)


(* ::Input:: *)
(*fit28["BestFitParameters"]*)
(*nd={{siInterp[300],0.0},{siO2Interp[300],16.0},{si3n4Interp[300],1090.0},{alInterp[300],130.0},{alF3Interp[300],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[300],0.0},{alInterp[300],130.0},{alF3Interp[300],85.0},{1.0,0.0}};*)
(*lam=300;*)
(*fractionS=0.9;*)
(*r30=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r30b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t30, r30d}=Transpose[r30];*)
(*{t30, r30db}=Transpose[r30b];*)
(*r30dif=Transpose[{t30,(r30d-r30db)}];*)
(*ListLogPlot[{r30,r30b}, Joined->True]*)
(*ListPlot[r30dif]*)
(*nd={{si3n4Interp[300],0.0},{alInterp[300],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d30=Drop[Import["S7W30.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d30];*)
(*weights=1/wt1^2;*)
(*fit30=NonlinearModelFit[d30,parrattR[nd,theta,lam,fractionS],{{th,103.832}, {ndx,0.937094}, {beta,0.06617}}, {theta},Weights->weights];*)
(*fit30["ParameterTable"]*)
(*Show[LogPlot[fit30[theta],{theta,0,80}],ListLogPlot[d30]]*)


(* ::Section::Closed:: *)
(*32*)


(* ::Input:: *)
(*compute[396,397,32,10]*)


(* ::Input:: *)
(*fit30["BestFitParameters"]*)
(*nd={{siInterp[320],0.0},{siO2Interp[320],16.0},{si3n4Interp[320],1090.0},{alInterp[320],130.0},{alF3Interp[320],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[320],0.0},{alInterp[320],130.0},{alF3Interp[320],85.0},{1.0,0.0}};*)
(*lam=320;*)
(*fractionS=0.9;*)
(*r32=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r32b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t32, r32d}=Transpose[r32];*)
(*{t32, r32db}=Transpose[r32b];*)
(*r32dif=Transpose[{t32,(r32d-r32db)}];*)
(*ListLogPlot[{r32,r32b}, Joined->True]*)
(*ListPlot[r32dif]*)
(*nd={{si3n4Interp[320],0.0},{alInterp[320],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d32=Drop[Import["S7W32.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d32];*)
(*weights=1/wt1^2;*)
(*fit32=NonlinearModelFit[d32,parrattR[nd,theta,lam,fractionS],{{th,102.683}, {ndx,0.931}, {beta,0.0776716}}, {theta},Weights->weights];*)
(*fit32["ParameterTable"]*)
(*Show[LogPlot[fit32[theta],{theta,0,80}],ListLogPlot[d32]]*)


(* ::Section::Closed:: *)
(*34*)


(* ::Input:: *)
(*compute[398,399,34,10]*)


(* ::Input:: *)
(*fit32["BestFitParameters"]*)
(*nd={{siInterp[340],0.0},{siO2Interp[340],16.0},{si3n4Interp[340],1090.0},{alInterp[340],130.0},{alF3Interp[340],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[340],0.0},{alInterp[340],130.0},{alF3Interp[340],85.0},{1.0,0.0}};*)
(*lam=340;*)
(*fractionS=0.9;*)
(*r34=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r34b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t34, r34d}=Transpose[r34];*)
(*{t34, r34db}=Transpose[r34b];*)
(*r34dif=Transpose[{t34,(r34d-r34db)}];*)
(*ListLogPlot[{r34,r34b}, Joined->True]*)
(*ListPlot[r34dif]*)
(*nd={{si3n4Interp[340],0.0},{alInterp[340],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d34=Drop[Import["S7W32.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d34];*)
(*weights=1/wt1^2;*)
(*fit34=NonlinearModelFit[d34,parrattR[nd,theta,lam,fractionS],{{th,98.3204}, {ndx,0.930223}, {beta,0.091346}}, {theta},Weights->weights];*)
(*fit34["ParameterTable"]*)
(*Show[LogPlot[fit34[theta],{theta,0,80}],ListLogPlot[d34]]*)


(* ::Section:: *)
(*36*)


(* ::Input:: *)
(*compute[400,401,36,10]*)


(* ::Input:: *)
(*fit34["BestFitParameters"]*)
(*nd={{siInterp[360],0.0},{siO2Interp[360],16.0},{si3n4Interp[360],1090.0},{alInterp[360],130.0},{alF3Interp[360],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[360],0.0},{alInterp[360],130.0},{alF3Interp[360],85.0},{1.0,0.0}};*)
(*lam=360;*)
(*fractionS=0.9;*)
(*r36=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r36b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t36, r36d}=Transpose[r36];*)
(*{t36, r36db}=Transpose[r36b];*)
(*r36dif=Transpose[{t36,(r36d-r36db)}];*)
(*ListLogPlot[{r36,r36b}, Joined->True]*)
(*ListPlot[r36dif]*)
(*nd={{si3n4Interp[360],0.0},{alInterp[360],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d36=Drop[Import["S7W36.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d36];*)
(*weights=1/wt1^2;*)
(*fit36=NonlinearModelFit[d36,parrattR[nd,theta,lam,fractionS],{{th,116.181}, {ndx,0.93198}, {beta,0.10576}}, {theta},Weights->weights];*)
(*fit36["ParameterTable"]*)
(*Show[LogPlot[fit36[theta],{theta,0,80}],ListLogPlot[d36]]*)


(* ::Section:: *)
(*38*)


(* ::Input:: *)
(*compute[402,403,38,10]*)
(*(*Look at the two data seperatly to solve for the problems in the middle*)*)


(* ::Input:: *)
(*fit36["BestFitParameters"]*)
(*nd={{siInterp[380],0.0},{siO2Interp[380],16.0},{si3n4Interp[380],1090.0},{alInterp[380],130.0},{alF3Interp[380],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[380],0.0},{alInterp[380],130.0},{alF3Interp[380],85.0},{1.0,0.0}};*)
(*lam=380;*)
(*fractionS=0.9;*)
(*r38=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r38b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t38, r38d}=Transpose[r38];*)
(*{t38, r38db}=Transpose[r38b];*)
(*r38dif=Transpose[{t38,(r38d-r38db)}];*)
(*ListLogPlot[{r38,r38b}, Joined->True]*)
(*ListPlot[r38dif]*)
(*nd={{si3n4Interp[380],0.0},{alInterp[380],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d38=Drop[Import["S7W38.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d38];*)
(*weights=1/wt1^2;*)
(*fit38=NonlinearModelFit[d38,parrattR[nd,theta,lam,fractionS],{{th,75.9811}, {ndx,0.934414}, {beta,0.129498}}, {theta},Weights->weights];*)
(*fit38["ParameterTable"]*)
(*Show[LogPlot[fit38[theta],{theta,0,80}],ListLogPlot[d38]]*)


(* ::Input:: *)
(*(*The funny points are throwing the fit off*)*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*41*)


(* ::Input:: *)
(*{\[Theta],iDatan1,c3,c4}=readRun[404];*)
(**)
(*{\[Theta]2,iDatan2,c3,c4}=readRun[405]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iDatan1~Join~(iDatan2*10^-(10-9)); *)
(**)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[41];*)
(*refla=(iData-d9)*10^(-9)/((i0-d9)*10^(-9));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[41];*)
(*Export["S7W"<>ToString[41]<>".txt",Transpose[{\[Theta],refla}],"Table"];*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refla}]]*)


(* ::Input:: *)
(*fit38["BestFitParameters"]*)
(*nd={{siInterp[410],0.0},{siO2Interp[410],16.0},{si3n4Interp[410],1090.0},{alInterp[410],130.0},{alF3Interp[410],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[410],0.0},{alInterp[410],130.0},{alF3Interp[410],85.0},{1.0,0.0}};*)
(*lam=410;*)
(*fractionS=0.9;*)
(*r41=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r41b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t41, r41d}=Transpose[r41];*)
(*{t41, r41db}=Transpose[r41b];*)
(*r41dif=Transpose[{t41,(r41d-r41db)}];*)
(*ListLogPlot[{r41,r41b}, Joined->True]*)
(*ListPlot[r41dif]*)
(*nd={{si3n4Interp[410],0.0},{alInterp[410],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d41=Drop[Import["S7W41.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d41];*)
(*weights=1/wt1^2;*)
(*fit41=NonlinearModelFit[d41,parrattR[nd,theta,lam,fractionS],{{th,88.3133}, {ndx,0.909345}, {beta,0.160618}}, {theta},Weights->weights];*)
(*fit41["ParameterTable"]*)
(*Show[LogPlot[fit41[theta],{theta,0,80}],ListLogPlot[d41]]*)


(* ::Section:: *)
(*44*)


(* ::Input:: *)
(*{\[Theta],iDatan1,c3,c4}=readRun[406];*)
(**)
(*{\[Theta]2,iDatan2,c3,c4}=readRun[407]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iDatan1~Join~(iDatan2*10^-(1)); *)
(*(*Check data book for correct gains, this is incorrect*)*)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[44];*)
(*refla=(iData-d9)*10^(-9)/((i0-d9)*10^(-9));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[44];*)
(*Export["S7W"<>ToString[44]<>".txt",Transpose[{\[Theta],refla}],"Table"];*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refla}]]*)


(* ::Input:: *)
(*fit41["BestFitParameters"]*)
(*nd={{siInterp[440],0.0},{siO2Interp[440],16.0},{si3n4Interp[440],1090.0},{alInterp[440],130.0},{alF3Interp[440],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[440],0.0},{alInterp[440],130.0},{alF3Interp[440],85.0},{1.0,0.0}};*)
(*lam=440;*)
(*fractionS=0.9;*)
(*r44=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r44b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t44, r44d}=Transpose[r44];*)
(*{t44, r44db}=Transpose[r44b];*)
(*r44dif=Transpose[{t44,(r44d-r44db)}];*)
(*ListLogPlot[{r44,r44b}, Joined->True]*)
(*ListPlot[r44dif]*)
(*nd={{si3n4Interp[440],0.0},{alInterp[440],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d44=Drop[Import["S7W44.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d44];*)
(*weights=1/wt1^2;*)
(*fit44=NonlinearModelFit[d44,parrattR[nd,theta,lam,fractionS],{{th,205.84}, {ndx,0.67547}, {beta,0.0205868}}, {theta},Weights->weights];*)
(*fit44["ParameterTable"]*)
(*Show[LogPlot[fit44[theta],{theta,0,80}],ListLogPlot[d44]]*)


(* ::Section:: *)
(*47*)


(* ::Input:: *)
(*{\[Theta],iDatan1,c3,c4}=readRun[408];*)
(**)
(*{\[Theta]2,iDatan2,c3,c4}=readRun[409]; *)
(*\[Theta]=\[Theta]1~Join~\[Theta]2;*)
(*iData=iDatan1~Join~(iDatan2*10^-(10-8)); *)
(**)
(*i0Interp=Interpolation[Transpose[{\[Lambda],inaught}]];*)
(*i0=i0Interp[47];*)
(*refla=(iData-d8)*10^(-8)/((i0-d8)*10^(-8));*)
(*linInterp=Interpolation[Transpose[{\[Lambda],inaught}],InterpolationOrder->1];*)
(*i0lin=linInterp[47];*)
(*Export["S7W"<>ToString[47]<>".txt",Transpose[{\[Theta],refla}],"Table"];*)
(**)
(*ListLogPlot[Transpose[{\[Theta],refla}]]*)


(* ::Input:: *)
(*fit44["BestFitParameters"]*)
(*nd={{siInterp[470],0.0},{siO2Interp[470],16.0},{si3n4Interp[470],1090.0},{alInterp[470],130.0},{alF3Interp[470],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[470],0.0},{alInterp[470],130.0},{alF3Interp[470],85.0},{1.0,0.0}};*)
(*lam=470;*)
(*fractionS=0.9;*)
(*r47=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r47b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t47, r47d}=Transpose[r47];*)
(*{t47, r47db}=Transpose[r47b];*)
(*r47dif=Transpose[{t47,(r47d-r47db)}];*)
(*ListLogPlot[{r47,r47b}, Joined->True]*)
(*ListPlot[r47dif]*)
(*nd={{si3n4Interp[470],0.0},{alInterp[470],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d47=Drop[Import["S7W47.txt","Table"],7];*)
(*{ang,wt1}=Transpose[d47];*)
(*weights=1/wt1^2;*)
(*fit47=NonlinearModelFit[d47,parrattR[nd,theta,lam,fractionS],{{th,238.762}, {ndx,0.65047}, {beta,0.0469166}}, {theta},Weights->weights];*)
(*fit47["ParameterTable"]*)
(*Show[LogPlot[fit47[theta],{theta,0,80}],ListLogPlot[d47]]*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Plotting Thicknesses and Constants*)


(* ::Input:: *)
(*th18=th/.fit18["BestFitParameters"][[1]];*)
(*th20=th/.fit20["BestFitParameters"][[1]];*)
(*th22=th/.fit22["BestFitParameters"][[1]];*)
(*th24=th/.fit24["BestFitParameters"][[1]];*)
(*th26=th/.fit26["BestFitParameters"][[1]];*)
(*th28=th/.fit28["BestFitParameters"][[1]];*)
(*th30=th/.fit30["BestFitParameters"][[1]];*)
(*th32=th/.fit32["BestFitParameters"][[1]];*)
(*th34=th/.fit34["BestFitParameters"][[1]];*)
(*th36=th/.fit36["BestFitParameters"][[1]];*)
(*th38=th/.fit38["BestFitParameters"][[1]];*)
(*th41=th/.fit41["BestFitParameters"][[1]];*)
(*th44=th/.fit44["BestFitParameters"][[1]];*)
(*th47=th/.fit47["BestFitParameters"][[1]];*)
(*thicknesses={th18,th20,th22,th24,th26,th28,th30,th32,th34,th36,th38,th41,th44,th47};*)
(*ListPlot[thicknesses]*)
(*Mean[thicknesses]*)


(* ::Input:: *)
(*n18=ndx/.fit18["BestFitParameters"][[2]];*)
(*n20=ndx/.fit20["BestFitParameters"][[2]];*)
(*n22=ndx/.fit22["BestFitParameters"][[2]];*)
(*n24=ndx/.fit24["BestFitParameters"][[2]];*)
(*n26=ndx/.fit26["BestFitParameters"][[2]];*)
(*n28=ndx/.fit28["BestFitParameters"][[2]];*)
(*n30=ndx/.fit30["BestFitParameters"][[2]];*)
(*n32=ndx/.fit32["BestFitParameters"][[2]];*)
(*n34=ndx/.fit34["BestFitParameters"][[2]];*)
(*n36=ndx/.fit36["BestFitParameters"][[2]];*)
(*n38=ndx/.fit38["BestFitParameters"][[2]];*)
(*n41=ndx/.fit41["BestFitParameters"][[2]];*)
(*n44=ndx/.fit44["BestFitParameters"][[2]];*)
(*n47=ndx/.fit47["BestFitParameters"][[2]];*)
(*ns={n18,n20,n22,n24,n26,n28,n30,n32,n34,n36,n38,n41,n44,n47};*)
(*(**)*)
(*ListPlot[ns]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*k18=beta/.fit18["BestFitParameters"][[3]];*)
(*k20=beta/.fit20["BestFitParameters"][[3]];*)
(*k22=beta/.fit22["BestFitParameters"][[3]];*)
(*k24=beta/.fit24["BestFitParameters"][[3]];*)
(*k26=beta/.fit26["BestFitParameters"][[3]];*)
(*k28=beta/.fit28["BestFitParameters"][[3]];*)
(*k30=beta/.fit30["BestFitParameters"][[3]];*)
(*k32=beta/.fit32["BestFitParameters"][[3]];*)
(*k34=beta/.fit34["BestFitParameters"][[3]];*)
(*k36=beta/.fit36["BestFitParameters"][[3]];*)
(*k38=beta/.fit38["BestFitParameters"][[3]];*)
(*k41=beta/.fit41["BestFitParameters"][[3]];*)
(*k44=beta/.fit44["BestFitParameters"][[3]];*)
(*k47=beta/.fit47["BestFitParameters"][[3]];*)
(*ks={k18,k20,k22,k24,k26,k28,k30,k32,k34,k36,k38,k41,k44,k47};*)
(*ListPlot[ks]*)
(*(*Go through peices to make it work*)*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Fitting AlF3*)


(* ::Section:: *)
(*18 nm*)


(* ::Input:: *)
(*nd={{siInterp[180],0.0},{siO2Interp[180],18.0},{si3n4Interp[180],1090.0},{alInterp[180],130.0},{alF3Interp[180],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[180],0.0},{alInterp[180],130.0},{alF3Interp[180],85.0},{1.0,0.0}};*)
(*lam=180;*)
(*fractionS=0.95;*)
(*r18=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r18b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t18, r18d}=Transpose[r18];*)
(*{t18, r18db}=Transpose[r18b];*)
(*r18dif=Transpose[{t18,(r18d-r18db)}];*)
(*ListLogPlot[{r18,r18b}, Joined->True]*)
(*alF3Interp[180];*)
(*(*Include SiO2 and Si in the nd fit. Si for substrate with SiO2 layer below, just like above.*)*)
(*nd={{si3n4Interp[180],1100.0},{alInterp[180],th},{}{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d18=Drop[Import["S7W18M.txt","Table"],7];*)
(*{ang,wt}=Transpose[d18];*)
(*weights=1/wt^2;*)
(*(*Need to calculate fractionS with new code. Between .85 and .95*)*)
(*(*fit18=NonlinearModelFit[d18,parrattR[nd,theta,lam,fractionS],{{th,112.0}, {ndx,0.951}, {beta,0.0205}}, {theta}];*)
(*fit18["ParameterTable"]*)*)
(*fit18=NonlinearModelFit[d18,parrattR[nd,theta,lam,fractionS],{{th,112.0}, {ndx,0.951}, {beta,0.0205}}, {theta},Weights->weights];*)
(*fit18["ParameterTable"]*)
(*Show[LogPlot[fit18[theta],{theta,0,80}],ListLogPlot[d18]]*)


(* ::Input:: *)
(**)


(* ::Section::Closed:: *)
(*20*)


(* ::Input:: *)
(*nd={{siInterp[200],0.0},{siO2Interp[200],16.0},{si3n4Interp[200],1090.0},{alInterp[200],130.0},{alF3Interp[200],85.0},{1.0,0.0}};*)
(*nd2={{si3n4Interp[200],0.0},{alInterp[200],130.0},{alF3Interp[200],85.0},{1.0,0.0}};*)
(*lam=200;*)
(*fractionS=0.9;*)
(*r20=Table[{theta,parrattR[nd,theta,lam,fractionS]},{theta,0,80,1}];*)
(*r20b=Table[{theta,parrattR[nd2,theta,lam,fractionS]},{theta,0,80,1}];*)
(*{t20, r20d}=Transpose[r20];*)
(*{t20, r20db}=Transpose[r20b];*)
(*r20dif=Transpose[{t20,(r20d-r20db)}];*)
(*ListLogPlot[{r20,r20b}, Joined->True]*)
(*ListPlot[r20dif]*)
(*alF3Interp[200];*)
(*nd={{si3n4Interp[200],1100.0},{alInterp[200],177.033},{alF3Interp[200],th},{ndx+beta*I,85.0},{1.0,0.0}};*)
(*d201=Take[Import["S7W20M.txt","Table"],{8,80}];*)
(*d202=Take[Import["S7W20M.txt","Table"],{86,160}];*)
(*d20=d201~Join~d202;*)
(*{ang,wt}=Transpose[d20];*)
(*weights=1/wt^2;*)
(*(*fit20=NonlinearModelFit[d20,parrattR[nd,theta,lam,fractionS],{{th,130.0}, {ndx,0.936}, {beta,0.102}}, {theta}];*)
(*fit20["ParameterTable"]*)*)
(*fit20=NonlinearModelFit[d20,parrattR[nd,theta,lam,fractionS],{{th,104.536}, {ndx,0.963}, {beta,0.036}}, {theta},Weights->weights];*)
(*fit20["ParameterTable"]*)
(*Show[LogPlot[fit20[theta],{theta,0,80}],ListLogPlot[d20]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Section:: *)
(**)


(* ::Input:: *)
(*ndx={9,10,11,12,11,9}*)
(*ListPlot[ndx]*)


(* ::Input:: *)
(*lam={21,23,25,27,29,31}*)
(*ndxlam=Transpose[{lam,ndx}]*)


(* ::Input:: *)
(*ListPlot[ndxlam]*)
