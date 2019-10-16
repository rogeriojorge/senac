# senac
Stellarator Equilibrium Near-Axis Code


Getting Started

The executable can be found in senac.sh, where all the input parameters are given. To run it, make it executable by writing chmod +x senac.sh and then run it using ./senac.sh
The file senac.sh is a bash script, first using Mathematica to read from a given input file, construct the near-axis surface of constant psi, and plotting it. Then it is able to run VMEC using the resulting surface and REGCOIL.


There are two different ways to run SENAC:

1 - Give a project_name in proj="project_name" different from a folder present in the VMEC folder, e.g., proj="test". The code will then read the parameters from the surf_input.txt file. If the mu variable is present, it reads mu, delta and B0 from surf_input. Otherwise, e.g., if mua is written in surf_input.txt instead of mu, it reads the Fourier coefficients in the file and then fits the parameters of the near-axis expansion to that surface.
2 - Give a project name in proj="project_name" equal to a VMEC output file that is in one of the folders in VMEC, e.g., proj="W7X". SENAC will then fit the parameters of the near-axis expansion to the innermost surface of that VMEC file.
 

Prerequisites

Mathematica
Add wolframscript to the workspace, or give its location to senac.sh to run SENAC
For first time Mac users, run 
sudo ln -s /Applications/Mathematica.app/Contents/MacOS/WolframScript /usr/local/bin/
sudo ln -s /Applications/MATLAB_RXXXx.app/bin/matlab /usr/local/bin

VMEC (optional)

https://princetonuniversity.github.io/STELLOPT/STELLOPT

REGCOIL (optional)

https://github.com/landreman/regcoil
