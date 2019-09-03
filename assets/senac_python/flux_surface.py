import numpy as np
import matplotlib.pyplot as plt
from input_senac import NFP, B0s, B0c, deltal, deltas, deltac, mul, mus, muc, nPointsPhiGamma, nPointsPhi, nPointsOmega
from frenetserret import tors, sprime, k0, t0, axis
from scipy import integrate

## Lowest Order Psi - rho^2
# Magnetic Field on the Axis
def B0(phi):
	B0temp = 0
	i = 0
	for i in range(len(B0s)):
		B0temp = B0temp + B0s[i]*np.sin(i*NFP*phi)
	j = 0
	for j in range(len(B0c)):
		B0temp = B0temp + B0c[j]*np.cos(j*NFP*phi)
	return B0temp

# Rotation Angle delta of the ellipse
def delta(phi):
	deltatemp = deltal*phi
	i = 0
	for i in range(len(deltas)):
		deltatemp = deltatemp - deltas[i]*np.sin(i*NFP*phi)
	j = 0
	for j in range(len(deltac)):
		deltatemp = deltatemp + deltac[j]*np.cos(j*NFP*phi)
	return deltatemp

# Modified Eccentricity mu of the ellipse
def mu(phi):
	mutemp = mul*phi
	i = 0
	for i in range(len(mus)):
		mutemp = mutemp - mus[i]*np.sin(i*NFP*phi)
	j = 0
	for j in range(len(muc)):
		mutemp = mutemp + muc[j]*np.cos(j*NFP*phi)
	return mutemp

# Integrated Torsion gamma of the axis
def gamma(phi):
	if phi < 2*np.pi/nPointsPhiGamma:
	#if phi < 4*np.pi:
		return 0
	else:
		phiTable = np.linspace(0., phi, num=np.floor(phi*nPointsPhiGamma/(2*np.pi)))
		integrand = np.zeros((len(phiTable)))
		for i in range(len(phiTable)):
			integrand[i] = tors(phiTable[i])*sprime(phiTable[i])
		gammaTemp = integrate.simps(integrand, phiTable)
		return gammaTemp

## Psi - rho^3
# To be implemented

## Calculate surface of constant psi
# Radius for each Flux Surface from Mercier's Theory
def rho(phi, omega, psi, ngamma):
	rhotemp = np.sqrt(psi/B0(phi))*((1-mu(phi)**2)**0.25)/np.sqrt(1+mu(phi)*np.cos(2*(omega-ngamma+delta(phi))))
	return rhotemp

# Position Vector of the Surface
def rPos(phi, omega, psi, ngamma):
    rPostemp = axis(phi)+rho(phi,omega,psi,ngamma)*(np.cos(omega-ngamma)*k0(phi)+np.sin(omega-ngamma)*t0(phi))
    return rPostemp

# Calculate 3D position vector r(phi, omega) from the definitions above
def nrPos(psi, maxtheta, maxphi):
	omegaTable = np.linspace(0., maxtheta, num=nPointsOmega)
	phiTable = np.linspace(0., maxphi, num=nPointsPhi)
	nrPosTemp = np.zeros((len(phiTable),len(omegaTable),3))
	for i in range(len(phiTable)):
		ngamma = gamma(phiTable[i])
		for j in range(len(omegaTable)):
			nrPosTemp[i][j] = np.array(rPos(phiTable[i], omegaTable[j], psi, ngamma))
	x = nrPosTemp[0:nPointsPhi+1,0:nPointsOmega+1,0]
	y = nrPosTemp[0:nPointsPhi+1,0:nPointsOmega+1,1]
	z = nrPosTemp[0:nPointsPhi+1,0:nPointsOmega+1,2]
	return x, y, z