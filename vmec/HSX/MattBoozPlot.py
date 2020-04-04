#!/usr/bin/env python

print "If any arguments are pdf, a pdf will be saved."

#filenames = ['20180404-01-015_R0c_0.045_Z0s_0.045_B1c_0.9_A10/boozmn_QS.nc', \
#                 '20180404-01-018_R0c_0.045_Z0s_0.045_B1c_0.9_A80/boozmn_QS.nc']

filenames = ['boozmn_wout_HSX.nc']
#filenames = ['boozmn_partialQA_r0.1_finiteRNonlinearAltSplines_nphi201_ntheta40_mpol12_ntor12_ns101.nc']

aspect_ratios = ['10','80']

max_m = 5
max_n = 5

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
import sys
import math
import inspect
import os
from os.path import abspath

#if len(sys.argv) != 2:
#    print "Error! You must specify 1 argument: the boozmnXXX.nc file."
#    exit(1)

makePDF = False
for arg in sys.argv:
   if arg.lower()=='pdf':
      makePDF = True

#fig = plt.figure(figsize=(9,4))
fig = plt.figure(figsize=(4.5,4))
fig.patch.set_facecolor('white')


for which_subplot in range(1):
    #plt.subplot(1,2,which_subplot+1)

    #filename = sys.argv[1]
    filename = filenames[which_subplot]
    print "About to try loading file ",filename
    f = netcdf.netcdf_file(filename,mode='r',mmap=False)

    phi_b = f.variables['phi_b'][()]
    ns_b = f.variables['ns_b'][()]
    nfp_b = f.variables['nfp_b'][()]
    ixn_b = f.variables['ixn_b'][()]
    ixm_b = f.variables['ixm_b'][()]
    bmnc_b = f.variables['bmnc_b'][()]
    jlist = f.variables['jlist'][()]
    f.close()
    nmodes = len(ixn_b)

    s = (jlist-1.5)/(ns_b-1.0)

    bmnc_b = -bmnc_b # Flip signs since the plot looks nicer with the QA mode positive.

    backgroundColor='b'
    #QAColor=[0,0.7,0]
    QAColor=[1,0,0]
    mirrorColor=[0.7,0.5,0]
    m1color = 'b'
    m2color = 'm'
    m3color = 'darkgreen'
    #helicalColor=[1,0,1]

    scale_factor = np.max(np.abs(bmnc_b))

    # First, plot just the 1st mode of each type, so the legend looks nice.

    #for imode in range(nmodes):
    #    if ixn_b[imode]==0 and ixm_b[imode]==0:
    #        plt.plot(s,abs(bmnc_b[:,imode])/scale_factor, color=backgroundColor,label='m = 0, n = 0 (Background)')
    #        break
    largest_amplitude = 0
    largest_amplitude_m = 999
    largest_amplitude_n = 999
    for imode in range(nmodes):
        if np.abs(bmnc_b[-1,imode]) > largest_amplitude:
            if ixm_b[imode] == 0 and ixn_b[imode] == 0:
                continue
            else:
                largest_amplitude = np.abs(bmnc_b[-1,imode])
                largest_amplitude_m = ixm_b[imode]
                largest_amplitude_n = ixn_b[imode]
    print "Largest symmetry-breaking mode at the edge has m=",largest_amplitude_m,", n=",largest_amplitude_n


    for imode in range(nmodes):
        if ixn_b[imode]==largest_amplitude_n and ixm_b[imode]==largest_amplitude_m:
            plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=m3color,label=r'Dominant Mode, $m = '+str(largest_amplitude_m)+'$, $n = '+str(largest_amplitude_n)+'$')
            break
    for imode in range(nmodes):
        if ixn_b[imode]==0 and ixm_b[imode]!=0:
                if ixn_b[imode]==largest_amplitude_n and ixm_b[imode]==largest_amplitude_n:
                    break
                plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=QAColor,label=r'Symmetry-breaking, $n = 0$')
                break
    for imode in range(nmodes):
        if ixn_b[imode]!=0 and ixm_b[imode]==0:
                if ixn_b[imode]==largest_amplitude_n and ixm_b[imode]==largest_amplitude_n:
                    break
                plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=mirrorColor,label=r'Symmetry-breaking, $m = 0$')
                break
    for imode in range(nmodes):
        if ixn_b[imode]!=0 and ixm_b[imode]==1:
                if ixn_b[imode]==largest_amplitude_n and ixm_b[imode]==largest_amplitude_n:
                    break
                plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=m1color,label=r'Symmetry-breaking, $m = 1$')
                break
    for imode in range(nmodes):
        if ixn_b[imode]!=0 and ixm_b[imode]>=2:
                if ixn_b[imode]==largest_amplitude_n and ixm_b[imode]==largest_amplitude_n:
                    break
                plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=m2color,label=r'Symmetry-breaking, $m \geq 2$')
                break
    #for imode in range(nmodes):
    #    if ixn_b[imode]!=0 and ixm_b[imode]==3:
    #        plt.plot(np.sqrt(s),(bmnc_b[:,imode])/scale_factor, color=m3color,label=r'Symmetry-breaking, m >= 3')
    #        break
    #plt.legend(fontsize='x-small',loc=2)
    #plt.legend(fontsize=9,loc=2,frameon=False)
    r_over_a_including_0 = np.insert(np.sqrt(s),0,0)
    eta_bar_over_aspect = 0.632 / 10
    predicted_color='salmon'
    #plt.plot(r_over_a_including_0, eta_bar_over_aspect*r_over_a_including_0,':',color=predicted_color,label='$B_{1,0}$ predicted by construction')
    
    plt.legend(fontsize=9,loc=2)
    # Now that the legend is node, plot all modes

    for imode in range(nmodes):
        value_at_0 = 0
        data = (bmnc_b[:,imode])/scale_factor
        if np.abs(ixm_b[imode]) > max_m:
            continue
        if np.abs(ixn_b[imode]) > max_n * nfp_b:
            continue
        if ixn_b[imode]==0:
            if ixm_b[imode]==0:
                mycolor = backgroundColor
                continue
            else:
                mycolor = QAColor
        else:
            if ixm_b[imode]==0:
                mycolor = mirrorColor
                value_at_0 = data[0]
            elif ixm_b[imode]==1:
                mycolor = m1color
                value_at_0 = 0
            else:
                mycolor = m2color
                value_at_0 = 0
            #elif ixm_b[imode]==2:
            #    mycolor = m2color
            #    value_at_0 = 0
        if ixm_b[imode]==largest_amplitude_m and ixn_b[imode]==largest_amplitude_n:
           mycolor = m3color
        
        plt.plot(r_over_a_including_0,np.insert(data,0,value_at_0), color=mycolor)

    #plt.plot(r_over_a_including_0, eta_bar_over_aspect*r_over_a_including_0,':',color=predicted_color,label='$B_{1,0}$ predicted by construction')
    #plt.xlabel('$r / a$ (Sqrt normalized toroidal flux)')
    plt.xlabel('$r / a = \sqrt{s}$')
    #plt.title('Aspect ratio '+aspect_ratios[which_subplot]+'\nFourier harmonics $B_{m,n}$ in Boozer coordinates')
    plt.title('Fourier harmonics $B_{m,n}$ in Boozer coordinates [T]    ')
    #plt.ylim([1e-4,1.2])
    plt.xlim([0,1])
    if (which_subplot==1):
        plt.ylim([-0.0015,0.012])

    #titleString = "Plot generated by "+ abspath(inspect.getfile(inspect.currentframe()))
    #plt.figtext(0.5,0.01,titleString,horizontalalignment='center',verticalalignment='bottom',fontsize=8)
    #plt.figtext(0.5,0.99,'File = '+filename,horizontalalignment='center',verticalalignment='top',fontsize=10)

#plt.subplots_adjust(top=0.88,bottom=0.11,left=0.06,right=0.98)
plt.subplots_adjust(top=0.94,bottom=0.11,left=0.1,right=0.97)

if makePDF:
    print "Saving PDF"
    plt.savefig(os.path.basename(__file__) + ".pdf")
else:
    plt.show()


