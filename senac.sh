#!/usr/bin/env bash
# version 0.1 - R. Jorge IREAP/UMD September 2019
proj="HSX";   # project name for input/output files, with vmec output VMEC_files/wout_"proj".nc
#======SENAC INPUT PARAMETERS=====
readVMEC=1;
vmecInput="/Users/rogeriojorge/Dropbox/senac/VMEC_files/VMEC_input_template.txt"; #template VMEC input file
vmecOutput="/Users/rogeriojorge/Dropbox/senac/VMEC_files/wout_${proj}.nc"; #VMEC output file to read
nthetaM=20;   #resolution in theta to compute Mercier angle
nphiM=35;     #resolution in phi to compute Mercier angle
deltac0=1.5;  #initial point for deltac0 betweeon 0 and pi
deltal0=-10.0;#initial point for deltal
deltalmin=0.0;#minimum deltal to help fit
deltalmax=0.0; #maximum deltal to help fit (put equal to deltalmin to leave -1.2*vmecNFP<deltal<1.2*vmecNFP)
muc0=-0.7;    #initial point for muc0
mucMin=-0.9   #minimum muc0 to help fit
mucMax=-0.3   #maximum muc0 to help fit
maxiterations=350; #max number of iterations during fit parameter
plotFit=1;    #Mathematica plots fit results
maxm=5;       #Maximum m to output to VMEC
maxn=7;       #Maximum n to output to VMEC
maxRecursTheta=30; #Theta resolution in numerical integration
maxRecursPhi=150;  #Phi resolution in numerical integration
#======VEMAC=====
runVMEC=1;
runVMECplotOriginal=1;
#======VEMAC Plot=====
runVMECplotFit=1;

#======START================
 echo "===================SENAC===================="
 echo "Stellarator Equilibrium Near-Axis Code"
 echo "============================================"
 #==========================

#======Print Run Parameters=====
#echo "Parameters: nthetaM $nthetaM nphiM $nphiM maxiterations $maxiterations maxm $maxm maxn $maxn maxRecursTheta $maxRecursTheta maxRecursPhi $maxRecursPhi"

#======RUN SENAC=====
	echo "-----------------------"
	echo "Running SENAC Mathematica"
	wolframscript -noprompt -script main.wls $readVMEC $vmecInput $vmecOutput $nthetaM $nphiM $deltac0 $deltal0 $deltalmin $deltalmax $muc0 $mucMin $mucMax $maxiterations $plotFit $maxm $maxn $maxRecursTheta $maxRecursPhi | tee data/senac_${proj}_output.txt
#======RUN VMEC=====
if (( $runVMEC == 1)); then
	echo "-----------------------"
	echo "Running VMEC from fit"
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
#======RUN VMECplot=====
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
#pkill -9 WolframScript