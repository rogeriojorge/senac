import time
import numpy as np
import matplotlib.pyplot as plt
from plot_options import coolfig, coolfig3D, coolfigSurf3D
from input_senac import nPointsPhi, psi0, out_axis, out_curvature, out_fSurface, out_sprime, out_torsion, maxtheta, maxphi
from frenetserret import curv, tors, axis, sprime
from flux_surface import nrPos

def output_command(name, function):
    tic = time.time()
    eval(function)
    toc = time.time()
    print('Outputed %s in %.2f seconds' % (name, toc-tic))

phiTable = np.linspace(0., 2*np.pi, num=nPointsPhi)
if out_axis == 1:
    output_command('axis', r"coolfig3D(axis(phiTable), '', '', '', r'magnetic_axis')")
if out_curvature == 1:
    output_command('curvature', r"coolfig(phiTable, curv, r'toroidal angle $\phi$', r'curvature $\kappa$', r'curvature')")
if out_fSurface == 1:
    output_command('fsurface', r"coolfigSurf3D(nrPos(psi0, maxtheta, maxphi), '', '', '', r'fSurface')")
if out_sprime == 1:
    output_command('sprime', r"coolfig(phiTable, sprime, r'toroidal angle $\phi$', r'arclength derivative $ds/d\phi$', r'sprime')")
if out_torsion == 1:
    output_command('torsion', r"coolfig(phiTable, tors, r'toroidal angle $\phi$', r'torsion $\tau$', r'torsion')")

#Output to VMEC

