import numpy as np
from scipy.misc import derivative
from input_senac import NFP, RAXIS, ZAXIS, dxfs

# Magnetic Axis Curve
def R0(phi):
	R0temp = 0
	i = 0
	for i in range(len(RAXIS)):
		R0temp = R0temp + RAXIS[i]*np.cos(i*NFP*phi)
	return R0temp
def Z0(phi):
	Z0temp = 0
	i = 0
	for i in range(len(ZAXIS)):
		Z0temp = Z0temp - ZAXIS[i]*np.sin(i*NFP*phi)
	return Z0temp
def axis(phi):
	return [R0(phi)*np.cos(phi),R0(phi)*np.sin(phi),Z0(phi)]

# Derivatives of the Axis
	# First Derivative
def daxis(phi):
	return [derivative(R0, phi, dx=dxfs)*np.cos(phi)-R0(phi)*np.sin(phi),\
			derivative(R0, phi, dx=dxfs)*np.sin(phi)+R0(phi)*np.cos(phi),\
			derivative(Z0, phi, dx=dxfs)]
	# Second Derivative
def ddaxis(phi):
	return [derivative(R0, phi, dx=dxfs, n=2, order=5)*np.cos(phi)-2*derivative(R0, phi, dx=dxfs)*np.sin(phi)-R0(phi)*np.cos(phi),\
			derivative(R0, phi, dx=dxfs, n=2, order=5)*np.sin(phi)+2*derivative(R0, phi, dx=dxfs)*np.cos(phi)-R0(phi)*np.sin(phi),\
			derivative(Z0, phi, dx=dxfs, n=2, order=5)]
	# Third Derivative
def dddaxis(phi):
	return [derivative(R0, phi, dx=dxfs, n=3, order=5)*np.cos(phi)-3*derivative(R0, phi, dx=dxfs, n=2, order=5)*np.sin(phi)-3*derivative(R0, phi, dx=dxfs)*np.cos(phi)+R0(phi)*np.sin(phi),\
			derivative(R0, phi, dx=dxfs, n=3, order=5)*np.sin(phi)+3*derivative(R0, phi, dx=dxfs, n=2, order=5)*np.cos(phi)-3*derivative(R0, phi, dx=dxfs)*np.sin(phi)-R0(phi)*np.cos(phi),\
			derivative(Z0, phi, dx=dxfs, n=3, order=5)]

# Curvature
def curv(phi):
	return (np.linalg.norm(np.cross(daxis(phi),ddaxis(phi))))/(np.linalg.norm(daxis(phi))**3)
# Torsion
def tors(phi):
	return (np.dot(np.cross(daxis(phi),ddaxis(phi)),dddaxis(phi)))/(np.linalg.norm(np.cross(daxis(phi),ddaxis(phi)))**2)

# ArcLength derivative for chain rule
def sprime(phi):
	return np.linalg.norm(daxis(phi))

# Frenet-Serret Unit Vectors
	# tangent
def b0(phi):
	return np.asarray(daxis(phi))/sprime(phi)
	# normal
def k0(phi):
	return (np.asarray(ddaxis(phi))*sprime(phi)-np.asarray(daxis(phi))*derivative(sprime, phi, dx=dxfs))/((sprime(phi)**3)*curv(phi))
	# binormal
def t0(phi):
	return np.cross(b0(phi),k0(phi))