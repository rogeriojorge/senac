## Main Input Parameters
NFP   =  3 # Number of field periods
RAXIS =  [1., 0.025] # Fourier Cosine components of the axis
ZAXIS =  [0., 0.025] # Fourier Sine components of the axis
psi0 = 0.2 # Flux Surface to compute

## Lowest Order Flux Surface Parameters
B0s = [0., 0.]     # Fourier Sine components of the magnetic field on axis
B0c = [1., 0.]     # Fourier Cosine components of the magnetic field on axis
deltal = 3.        # Linear part of delta proportional to the toroidal angle
deltas = [0., 0.0] # Fourier Sine components of delta
deltac = [0., 0.]  # Fourier Cosine components of delta
mul = 0.           # Linear part of mu proportional to the toroidal angle
mus = [0.0, 0.]    # Fourier Sine components of mu
muc = [0.5, 0.0]   # Fourier Cosine components of mu

## Higher Order Flux Surface Parameters
# To be implemented

## Numerical Differentiation/Integration Parameters
dxfs  =  1e-4        # Derivative step-size for Frenet-Serret Calculation
nPointsPhiGamma = 20 # Number of points in the integral when computing integrated torsion
nPointsPhi = 50      # Resolution in the toroidal angle
nPointsOmega = 50    # Resolution in the poloidal angle

## Quantities to output
out_curvature = 1 # Output Curvature Figure
out_torsion = 1   # Output Torsion Figure
out_sprime = 1    # Output Sprime Figure
out_axis = 1      # Output Magnetic Axis Figure
out_fSurface = 1  # Output Flux Surface Figure

# Output parameters
maxtheta = 6.28
maxphi   = 6.28