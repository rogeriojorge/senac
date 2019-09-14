#!/usr/bin/env bash
# version 0.1 - R. Jorge IREAP/UMD September 2019
proj="LHD"; # project name for input/output files, with vmec output vmec/wout_"proj".nc
#================
currentDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
surfInput=${currentDIR}"/surf_input.txt"; #input file with surface parameters
vmecInput=${currentDIR}"/vmec/vmec_input_template.txt"; #template VMEC input file
vmecOutput=${currentDIR}"/vmec/${proj}/wout_${proj}.nc"; #VMEC output file to read
#======SENAC=====
runSENAC=1;              #1-> runs SENAC mathematica
outputToVMEC=0;          #compute Fourier Modes and output to VMEC
plotFit=1;               #Mathematica plots fit results
plotOriginal=1;          #Mathematica plots original surface
#======VMEC=====
runVMECofFit=0;          #run VMEC for fit
#======REGCOIL=====
runREGCOILfit=0;         #run REGCOIL for fit
plotRegcoilFit=0;        #Mathematica plots coils for fit
runREGCOILoriginal=0;    #run REGCOIL for original/VMEC file
plotRegcoilOriginal=0;   #Mathematica plots coils for original surface
#======SENAC INPUT PARAMETERS=====
ordern=2;                #Near-Axis Expansion Order (has to be greater than 2)
nModes=0;                #number of fourier components in mu, delta and B0
nsurfaces=1;             #number of surfaces to read and compare from VMEC
nthetaM=15;              #resolution in theta for fit and Mercier's coordinates
nphiM=20;                #resolution in phi for fit and Mercier's coordinates
maxiterations=2000;      #max number of iterations for fit
deltac0=1.1;             #initial point for deltac0 betweeon -pi and pi
deltal0=1.0;             #initial point for deltal
deltalmin=0.0;           #minimum deltal to help fit
deltalmax=0.0;           #maximum deltal to help fit (put equal to deltalmin to leave -1.2*vmecNFP<deltal<1.2*vmecNFP)
muc0=0.4;                #initial point for muc0
mucMin=-0.9;             #minimum muc0 to help fit
mucMax=0.9;              #maximum muc0 to help fit
maxm=6;                  #Maximum poloidal Fourier mode m to output to VMEC
maxn=7;                  #Maximum toroidal Fourier mode n to output to VMEC
maxRecursTheta=30;       #Theta resolution in numerical integration in Mercier to VMEC
maxRecursPhi=300;        #Phi resolution in numerical integration
#======PLOTTING PARAMETERS=====
exportBFieldSurface=1;   #0 -> Don't export figure of magnetic field on surface, 1 -> Do
nPlotTheta=80;           #number of interpolating points in theta
nPlotPhi=120;            #number of interpolating points in phi
plotPointsFig=50;        #plotpoints for 3D figure
maxRecursPlot=2;         #max recursion for 3D figure
ImageSizePlot=700;       #image size for 3D figure
ImageResolutionPlot=300; #resolution for 3D figure
nfigsSurf=3;             #number of surfaces to plot in 3D figure
nPlots=4;                #number of poloidal plots to save
npointsPolPlots=35;      #number of points for poloidal plots
nthetapointsBsurface=30; #plot points in theta for magnetic field on surface
nphipointsBsurface=30;   #plot points in phi for magnetic field on surface
npointsContourPlotREGCOIL=60;   #number of points in contourplot when finding coil contours in REGCOIL
npointsInterpCoilPosREGCOIL=60; #number of points for theta grid in REGCOIL
interpOrderCoilPosREGCOIL=3;    #interpolation order for theta grid in REGCOIL
coilsPerHalfPeriod=3;    #number of coils per half period to plot
numHalfPeriodsToPlot=0;  #0 -> plots the whole stellarator
#=====REGCOIL INPUT PARAMTERS========
REGCOILtargetvalue=0.3;
REGCOILseparation=0.08;
#=====TO BE IMPLEMENTED=============
readFit=0;    		     #1 -> reads fit parameters from text file
quasisymmetry=0;         #1 -> run quasisymmetric SENAC

#======START================
 echo "===================SENAC===================="
 echo "Stellarator Equilibrium Near-Axis Code"
 echo "============================================"
 echo ${proj}" Stellarator"
 mkdir -p data/${proj}/Figures
#==========================

#======Print Run Parameters=====
#echo "Parameters: nthetaM $nthetaM nphiM $nphiM maxiterations $maxiterations maxm $maxm maxn $maxn maxRecursTheta $maxRecursTheta maxRecursPhi $maxRecursPhi"

#======TEST INPUT PARAMETERS=====
if (($ordern < 2)); then
	echo "The order of the axis expansion should be above two"
	exit 1;
fi

#======RUN SENAC=====
if (( $runSENAC == 1)); then
	echo "-----------------------"
	echo "Running SENAC Mathematica"
	rm -f data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
	wolframscript -noprompt -script main.wls $proj $surfInput $readFit $outputToVMEC $vmecInput $vmecOutput $ordern $nsurfaces $nthetaM $nphiM $deltac0 $deltal0 $deltalmin $deltalmax $muc0 $mucMin $mucMax $nModes $maxiterations $plotFit $plotOriginal $maxm $maxn $maxRecursTheta $maxRecursPhi $nPlotTheta $nPlotPhi $plotPointsFig $maxRecursPlot $ImageSizePlot $ImageResolutionPlot $nfigsSurf $nPlots $nthetapointsBsurface $nphipointsBsurface $npointsPolPlots $exportBFieldSurface | tee data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
fi
#======RUN VMEC=====
if (( $runVMECofFit == 1)); then
	echo "-----------------------"
	echo "Running VMEC from fit"
	./vmec/xvmec2000 data/${proj}/input.${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt | tee -a data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
	if test -f "threed1.${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt"; then
    	rm threed1.${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt; rm timings.txt; rm mercier.${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt; rm parvmecinfo.txt; rm wout_${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt.txt; rm jxbout_${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt.nc; rm fort.8
	fi
	if test -f "wout_${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt.nc"; then
		mv wout_${proj}_SENACtoVMEC_ordern${ordern}_nmodes${nModes}.txt.nc data/${proj}/wout_${proj}_senac_ordern${ordern}_nmodes${nModes}.nc
	fi
fi
#======RUN REGCOIL=====
if (( $runREGCOILoriginal == 1)); then
	echo "-----------------------"
	echo "Running REGCOIL from original VMEC file"
	cd regcoil
	mv regcoil_in.senac regcoil_in.senac_temp
	head -n 19 regcoil_in.senac_temp > regcoil_in.senac
	rm regcoil_in.senac_temp
	echo '  target_value = '$REGCOILtargetvalue >> regcoil_in.senac
	echo '  separation = '$REGCOILseparation >> regcoil_in.senac
	echo '  wout_filename = "'${vmecOutput}'"' >> regcoil_in.senac
	echo '/' >> regcoil_in.senac
	./regcoil regcoil_in.senac
	mv regcoil_out.senac.nc ../vmec/${proj}/regcoil_out_${proj}.nc
	rm -f nescin.out
	cd ..
fi
if (( $runREGCOILfit == 1)); then
	echo "-----------------------"
	echo "Running REGCOIL from fit"
	get_abs_filename() {
  		echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
	}
	woutFile=$(get_abs_filename "data/${proj}/wout_${proj}_senac_ordern${ordern}_nmodes${nModes}.nc")
	cd regcoil
	mv regcoil_in.senac regcoil_in.senac_temp
	head -n 19 regcoil_in.senac_temp > regcoil_in.senac
	rm regcoil_in.senac_temp
	echo '  target_value = '$REGCOILtargetvalue >> regcoil_in.senac
	echo '  separation = '$REGCOILseparation >> regcoil_in.senac
	echo '  wout_filename = "'${woutFile}'"' >> regcoil_in.senac
	echo '/' >> regcoil_in.senac
	./regcoil regcoil_in.senac | tee -a ../data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
	mv regcoil_out.senac.nc ../data/${proj}/regcoil_out_${proj}_senac_ordern${ordern}_nmodes${nModes}.nc
	rm -f nescin.out
	cd ..
fi
#======RUN REGCOIL PLOT=====
if ( [[ $plotRegcoilOriginal == 1 ]] ); then
	echo "-----------------------"
	echo "Running SENAC Plot REGCOIL for Original Surface"
	plotRegcoilSENAC=1;
	wolframscript -noprompt -script src/plot_regcoil.wls $proj $ordern $nModes $ImageSizePlot $ImageResolutionPlot $coilsPerHalfPeriod $numHalfPeriodsToPlot $plotRegcoilSENAC $npointsContourPlotREGCOIL $npointsInterpCoilPosREGCOIL $interpOrderCoilPosREGCOIL $plotPointsFig $maxRecursPlot | tee -a data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
fi

if ( [[ $plotRegcoilFit == 1 ]]  ); then
	echo "-----------------------"
	echo "Running SENAC Plot REGCOIL for Fit Surface"
	plotRegcoilSENAC=2;
	wolframscript -noprompt -script src/plot_regcoil.wls $proj $ordern $nModes  $ImageSizePlot $ImageResolutionPlot $coilsPerHalfPeriod $numHalfPeriodsToPlot $plotRegcoilSENAC $npointsContourPlotREGCOIL $npointsInterpCoilPosREGCOIL $interpOrderCoilPosREGCOIL $plotPointsFig $maxRecursPlot | tee -a data/${proj}/senac_${proj}_output_order${ordern}_nmodes${nModes}.txt
fi

#pkill -9 WolframKernel; pkill -9 WolframScript
