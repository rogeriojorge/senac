#!/usr/bin/env bash
# version 0.1 - R. Jorge IREAP/UMD September 2019
proj="ARIES";   # project name for input/output files, with vmec output vmec/wout_"proj".nc
#======SENAC=====
runSENAC=0;   #put to 0 if only for plotting
readVMEC=1;   #put to 1 if reading axis and surface from VMEC
outputToVMEC=1; #compute Fourier Modes and output to VMEC
#======VMEC=====
runVMECofFit=0;
#======REGCOIL=====
runREGCOILoriginal=0;
runREGCOILfit=0;
#======VMECplot====
VMECplotOriginal=1;
VMECplotFit=0;
REGCOILplotOriginal=0;
REGCOILplotFit=1;
#======SENAC INPUT PARAMETERS=====
vmecInput="/Users/rogeriojorge/Dropbox/senac/vmec/vmec_input_template.txt"; #template VMEC input file
vmecOutput="/Users/rogeriojorge/Dropbox/senac/vmec/${proj}/wout_${proj}.nc"; #VMEC output file to read
nthetaM=20;   #resolution in theta to compute Mercier angle
nphiM=35;     #resolution in phi to compute Mercier angle
deltac0=1.5;  #initial point for deltac0 betweeon 0 and pi
deltal0=-1.0;#initial point for deltal
deltalmin=0.0;#minimum deltal to help fit
deltalmax=0.0; #maximum deltal to help fit (put equal to deltalmin to leave -1.2*vmecNFP<deltal<1.2*vmecNFP)
muc0=0.5;     #initial point for muc0
mucMin=0.1;   #minimum muc0 to help fit
mucMax=0.9;   #maximum muc0 to help fit
nModes=0;     #number of fourier components in mu, delta and B0
maxiterations=550; #max number of iterations during fit parameter
plotFit=1;    #Mathematica plots fit results
maxm=5;       #Maximum m to output to VMEC
maxn=7;       #Maximum n to output to VMEC
maxRecursTheta=30; #Theta resolution in numerical integration
maxRecursPhi=150;  #Phi resolution in numerical integration
#======VMECplot INPUT PARAMETERS======
nplotTheta=80;
nplotThetaSurf=80;
nplotPhiSurf=220;
nplotNthetaBSurf=40;
nplotNphiBSurf=50;
coilsPerHalfPeriod=3;
thetaShift=0;
#=====REGCOIL INPUT PARAMTERS========
REGCOILseparation = 0.07

#======START================
 echo "===================SENAC===================="
 echo "Stellarator Equilibrium Near-Axis Code"
 echo "============================================"
 echo ${proj}" Stellarator"
 mkdir -p data/${proj}/Figures
#==========================

#======Print Run Parameters=====
#echo "Parameters: nthetaM $nthetaM nphiM $nphiM maxiterations $maxiterations maxm $maxm maxn $maxn maxRecursTheta $maxRecursTheta maxRecursPhi $maxRecursPhi"

#======RUN SENAC=====
if (( $runSENAC == 1)); then
	echo "-----------------------"
	echo "Running SENAC Mathematica"
	rm -f data/${proj}/senac_${proj}_output.txt
	wolframscript -noprompt -script main.wls $proj $readVMEC $vmecInput $vmecOutput $nthetaM $nphiM $deltac0 $deltal0 $deltalmin $deltalmax $muc0 $mucMin $mucMax $nModes $maxiterations $plotFit $outputToVMEC $maxm $maxn $maxRecursTheta $maxRecursPhi | tee data/${proj}/senac_${proj}_output.txt
fi
#======RUN VMEC=====
if (( $runVMECofFit == 1)); then
	echo "-----------------------"
	echo "Running VMEC from fit"
	./vmec/xvmec2000 data/${proj}/${proj}_SENACtoVMEC_input.txt | tee -a data/${proj}/senac_${proj}_output.txt
	if test -f "jxbout_txt.nc"; then
    	rm *.txt
		rm jxbout_txt.nc
		rm fort.8
	fi
	if test -f "wout_txt.nc"; then
		mv wout_txt.nc data/${proj}/wout_${proj}_senac.nc
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
	echo 'separation = ${REGCOILseparation}'
	echo 'wout_filename = "'${vmecOutput}'"' >> regcoil_in.senac
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
	woutFile=$(get_abs_filename "data/${proj}/wout_${proj}_senac.nc")
	cd regcoil
	mv regcoil_in.senac regcoil_in.senac_temp
	head -n 19 regcoil_in.senac_temp > regcoil_in.senac
	rm regcoil_in.senac_temp
	echo 'separation = ${REGCOILseparation}'
	echo 'wout_filename = "'${woutFile}'"' >> regcoil_in.senac
	echo '/' >> regcoil_in.senac
	./regcoil regcoil_in.senac | tee -a ../data/${proj}/senac_${proj}_output.txt
	mv regcoil_out.senac.nc ../data/${proj}/regcoil_out_${proj}_senac.nc
	rm -f nescin.out
	cd ..
fi
#======RUN VMECplot=====
if (( $VMECplotOriginal == 1)); then
	echo "-----------------------"
	echo "Running VMECPLOT Original"
 	mkdir -p vmec/${proj}/Figures
	./vmecPlot vmec/${proj}/wout_${proj}.nc $proj $nplotTheta $nplotPhiSurf $nplotThetaSurf $nplotNthetaBSurf $nplotNphiBSurf $REGCOILplotOriginal "vmec/${proj}/regcoil_out_${proj}.nc" $coilsPerHalfPeriod $thetaShift
fi
if (( $VMECplotFit == 1)); then
	echo "-----------------------"
	echo "Running VMECPlot Fit"
	./vmecPlot data/${proj}/wout_${proj}_senac.nc $proj $nplotTheta $nplotPhiSurf $nplotThetaSurf $nplotNthetaBSurf $nplotNphiBSurf $REGCOILplotFit "data/${proj}/regcoil_out_${proj}_senac.nc" $coilsPerHalfPeriod $thetaShift | tee -a data/${proj}/senac_${proj}_output.txt
fi

#pkill -9 WolframKernel
#pkill -9 WolframScript
