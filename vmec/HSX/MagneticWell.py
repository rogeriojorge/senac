#!/usr/bin/env python3

import numpy as np
from scipy.io import netcdf

directory = ''
filenames = ['wout_HSX.nc']

N = len(filenames)

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(14,7))

numRows=2
numCols=3
plotNum = 1

d2_volume_d_psi2s = []

for j in range(N):
    plt.subplot(numRows,numCols,j+1)

    filename = directory + filenames[j]
    f = netcdf.netcdf_file(filename,'r',mmap=False)
    ns = f.variables['ns'][()]
    phi = f.variables['phi'][()]
    gmnc = f.variables['gmnc'][()]
    f.close()

    psi_a = phi[-1] / (2*np.pi)

    print('gmnc.shape:',gmnc.shape)
    s_full = np.linspace(0,1,ns,endpoint=True)
    ds = s_full[1] - s_full[0]
    s_half = s_full[1:] - (ds/2.0)
    data = gmnc[1:,0]
    plt.plot(s_half, data, '.-')
    plt.xlabel('s')
    plt.ylabel(r'$gmnc_{0,0}$')

    degree = 7
    max_index=50
    p = np.polyfit(s_half[:max_index], gmnc[1:max_index+1,0], degree)
    print("p: ",p)
    plt.plot(s_half, np.polyval(p,s_half),'-')
    print(p[degree-1])

    # In the next line, the minus sign is from the fact that sign(Jacobian) is <0.
    d2_volume_d_psi2 = -1/(psi_a*psi_a)*p[degree-1]
    print("Magnetic well, for plotting:",d2_volume_d_psi2)
    d2_volume_d_psi2s.append(d2_volume_d_psi2)

print("Final array of d2_volume_d_psi2:",d2_volume_d_psi2s)

plt.tight_layout()

plt.show()
