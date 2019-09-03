#!/usr/bin/env bash
# version 0.1 - R. Jorge IREAP/UMD September 2019
proj="LHD";   # project name for input/output files, with vmec output vmec/wout_"proj".nc
#======SENAC INPUT PARAMETERS=====
runSENAC=0;   #put to 0 if only for plotting
readVMEC=1;
vmecInput="/Users/rogeriojorge/Dropbox/senac/vmec/vmec_input_template.txt"; #template VMEC input file
vmecOutput="/Users/rogeriojorge/Dropbox/senac/vmec/wout_${proj}.nc"; #VMEC output file to read
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
plotFit=0;    #Mathematica plots fit results
maxm=4;       #Maximum m to output to VMEC
maxn=5;       #Maximum n to output to VMEC
maxRecursTheta=20; #Theta resolution in numerical integration
maxRecursPhi=100;  #Phi resolution in numerical integration
#======VMEC=====
runVMEC=0;
VMECplotOriginal=0;
VMECplotFit=1;
nplotTheta=10;
nplotPhiSurf=20;
nplotThetaSurf=20;
#======REGCOIL=====
runREGCOILoriginal=0;
runREGCOILfit=1;
REGCOILoriginalPlot=1;
REGCOILfitPlot=1;
coilsPerHalfPeriod=2;
thetaShift=0;

#======START================
 echo "===================SENAC===================="
 echo "Stellarator Equilibrium Near-Axis Code"
 echo "============================================"
 echo "Project name: "${proj}
 #==========================

#======Print Run Parameters=====
#echo "Parameters: nthetaM $nthetaM nphiM $nphiM maxiterations $maxiterations maxm $maxm maxn $maxn maxRecursTheta $maxRecursTheta maxRecursPhi $maxRecursPhi"

#======RUN SENAC=====
if (( $runSENAC == 1)); then
	echo "-----------------------"
	echo "Running SENAC Mathematica"
	wolframscript -noprompt -script main.wls $readVMEC $vmecInput $vmecOutput $nthetaM $nphiM $deltac0 $deltal0 $deltalmin $deltalmax $muc0 $mucMin $mucMax $maxiterations $plotFit $maxm $maxn $maxRecursTheta $maxRecursPhi | tee data/senac_${proj}_output.txt
fi
#======RUN VMEC=====
if (( $runVMEC == 1)); then
	echo "-----------------------"
	echo "Running VMEC from fit"
	./vmec/xvmec2000 data/wout_${proj}.nc_input.txt | tee -a data/senac_${proj}_output.txt
	if test -f "jxbout_txt.nc"; then
    	rm *.txt
		rm jxbout_txt.nc
		rm fort.8
	fi
	if test -f "wout_txt.nc"; then
		mv wout_txt.nc data/wout_${proj}_senac.nc | tee -a data/senac_${proj}_output.txt
	fi
fi
#======RUN REGCOIL=====
if (( $runREGCOILoriginal == 1)); then
	echo "-----------------------"
	echo "Running REGCOIL from original VMEC file"
	cd regcoil
	mv regcoil_in.senac regcoil_in.senac_temp
	head -n 21 regcoil_in.senac_temp > regcoil_in.senac
	rm regcoil_in.senac_temp
	echo 'wout_filename = "'${vmecOutput}'"' >> regcoil_in.senac
	echo '/' >> regcoil_in.senac
	./regcoil regcoil_in.senac
	mv regcoil_out.senac.nc regcoil_out_${proj}.nc
	cd ..
fi
if (( $runREGCOILfit == 1)); then
	echo "-----------------------"
	echo "Running REGCOIL from fit"
	get_abs_filename() {
  		# $1 : relative filename
  		echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
	}
	woutFile=$(get_abs_filename "data/wout_${proj}_senac.nc")
	cd regcoil
	mv regcoil_in.senac regcoil_in.senac_temp
	head -n 21 regcoil_in.senac_temp > regcoil_in.senac
	rm regcoil_in.senac_temp
	echo 'wout_filename = "'${woutFile}'"' >> regcoil_in.senac
	echo '/' >> regcoil_in.senac
	./regcoil regcoil_in.senac | tee -a ../data/senac_${proj}_output.txt
	mv regcoil_out.senac.nc ../data/regcoil_out_${proj}_senac.nc
	cd ..
fi
#======RUN VMECplot=====
if (( $VMECplotOriginal == 1)); then
	echo "-----------------------"
	echo "Running VMECPLOT Original"
	./vmecPlot vmec/wout_${proj}.nc $nplotTheta $nplotPhiSurf $nplotThetaSurf $REGCOILoriginalPlot "regcoil/regcoil_out_${proj}.nc" $coilsPerHalfPeriod $thetaShift
fi
if (( $VMECplotFit == 1)); then
	echo "-----------------------"
	echo "Running VMECPlot Fit"
	./vmecPlot data/wout_${proj}_senac.nc $nplotTheta $nplotPhiSurf $nplotThetaSurf $REGCOILfitPlot "data/regcoil_out_${proj}_senac.nc" $coilsPerHalfPeriod $thetaShift | tee -a data/senac_${proj}_output.txt
fi

#pkill -9 WolframKernel
#pkill -9 WolframScript
