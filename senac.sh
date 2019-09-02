#!/usr/bin/env bash
# version 0.1 - R. Jorge IREAP/UMD September 2019
#======START================
 echo "===================SENAC===================="
 echo "Stellarator Equilibrium Near-Axis Code"
 echo "============================================"
#======SENAC INPUT PARAMETERS=====
readVMEC=1;
vmecInput="/Users/rogeriojorge/Dropbox/senac/VMEC_files/VMEC_input_template.txt"; #template VMEC input file
proj="LHD";   # project name for file output from vmec output VMEC_files/wout_"proj".nc
vmecOutput="/Users/rogeriojorge/Dropbox/senac/VMEC_files/wout_${proj}.nc"; #VMEC output file to read
nthetaM=15;   #resolution in theta to compute Mercier angle
nphiM=25;     #resolution in phi to compute Mercier angle
deltac0=1.5; #initial point for deltac0
deltal0=-2.11; #initial point for deltal
muc0=0.5;    #initial point for muc0
mucMin=0.3   #minimum muc0 to help fit
maxiterations=100; #max number of iterations during fit parameter
perturbationscale=5.0 #Fit parameter
plotFit=1;    #Mathematica plots fit results
maxm=5;       #Maximum m to output to VMEC
maxn=6;       #Maximum n to output to VMEC
maxRecursTheta=30; #Theta resolution in numerical integration
maxRecursPhi=150;  #Phi resolution in numerical integration
#======VEMAC=====
runVMEC=1;
runVMECplotOriginal=0;
#======VEMAC Plot=====
runVMECplotFit=1;

#======Print Run Parameters=====
#echo "Parameters: nthetaM $nthetaM nphiM $nphiM maxiterations $maxiterations maxm $maxm maxn $maxn maxRecursTheta $maxRecursTheta maxRecursPhi $maxRecursPhi"

#======RUN SENAC=====
	echo "-----------------------"
	echo "Running SENAC Mathematica"
	wolframscript -noprompt -script main.wls $readVMEC $vmecInput $vmecOutput $nthetaM $nphiM $deltac0 $deltal0 $muc0 $mucMin $maxiterations $perturbationscale $plotFit $maxm $maxn $maxRecursTheta $maxRecursPhi
#======RUN VEMAC=====
if (( $runVMEC == 1)); then
	echo "-----------------------"
	echo "Running VMEC"
	./xvmec2000 data/wout_${proj}.nc_input.txt
	if test -f "jxbout_txt.nc"; then
    	rm *.txt
		rm jxbout_txt.nc
		rm fort.8
	fi
	if test -f "wout_txt.nc"; then
		mv wout_txt.nc data/wout_${proj}_senac.nc
	fi
fi
#======RUN VEMECplot=====
if (( $runVMECplotOriginal == 1)); then
	echo "-----------------------"
	echo "Running VMECPLOT Original"
	./vmecPlot VMEC_files/wout_${proj}.nc
fi
if (( $runVMECplotFit == 1)); then
	echo "-----------------------"
	echo "Running VMECPlot Fit"
	./vmecPlot data/wout_${proj}_senac.nc
fi

#pkill -9 WolframKernel