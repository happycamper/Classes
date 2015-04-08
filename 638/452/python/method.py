#!/usr/bin/python

import numpy as np
import scipy as sp
from scipy.integrate import quad
from scipy import linalg
import math
import cmath
import matplotlib.pyplot as plt


##Variables
subdivisions = 11
middle = (subdivisions + 1)/2
frequency = 300000000
SPEED_OF_LIGHT = 300000000
######


wavelength = SPEED_OF_LIGHT/frequency
k = 2.0*math.pi/wavelength
j = complex(0.0,1.0)
antenna_length = 0.47*wavelength
a = 0.005*wavelength
w = 2*math.pi*frequency
perm = 8.854e-12
r = .4418*wavelength
eta = 377

integrand_C = 1/(j*w*perm*4*math.pi)



deltaZ = antenna_length/(subdivisions) #in m



print "Voltage " + str(1/deltaZ)



V = np.zeros(shape=(subdivisions,1)) #voltage vector for excitation
V = np.complex128(V)
I = np.zeros(shape=(subdivisions,1))
I = np.complex128(I)
E = np.zeros(shape=(subdivisions,1))
E = np.complex128(E)
Ephi = np.zeros(shape=(subdivisions,1))
Ephi = np.complex128(E)
Z = np.zeros(shape=(subdivisions,subdivisions)) #zero out our impedance vector
Z = np.complex128(Z)
Zs = [] #a precalculated vector of distances along the antenna for quickness of computation
Angles = []

for F in range(0,240):
	Angles.append((math.pi*F)/120.0)

for x in range(1,subdivisions+1):
	if x == middle:
		V[middle-1,0] = complex((1/deltaZ),0.0)
	else:
		V[x-1,0] = complex(0.0,0.0)

	Zs.append((-1.0*antenna_length/2.0)+(deltaZ*(x-1))+(deltaZ/2.0))


def real_integrand(z,zm):
	R = math.sqrt(a**2+pow((zm-z),2))
	return sp.real((cmath.exp(-1.0*j*k*R)/pow(R,5))*((1+j*k*R)*((2*pow(R,2))-(3*a**3))+pow(k*a*R,2))*deltaZ)

def img_integrand(z,zm):
	R = math.sqrt(a**2+pow((zm-z),2))
	return sp.imag((cmath.exp(-1.0*j*k*R)/pow(R,5))*((1+j*k*R)*((2*pow(R,2))-(3*a**3))+pow(k*a*R,2))*deltaZ)

def real_eq2(z,zm):
	R = math.sqrt(a**2+pow((zm-z),2))
	return sp.real(-1.0*j*73.0*(deltaZ/wavelength)*((1+j*2*math.pi*R/wavelength)*(2-3*pow(a/R,2))+pow(2*math.pi*a/wavelength,2))*(math.cos(2*math.pi*R/wavelength)-j*math.sin(2*math.pi*R/wavelength))/(8*pow(math.pi,2)*wavelength*pow(R/wavelength,3)))

def imag_eq2(z,zm):
	R = math.sqrt(a**2+pow((zm-z),2))
	return sp.imag(-1.0*j*73.0*(deltaZ/wavelength)*((1+j*2*math.pi*R/wavelength)*(2-3*pow(a/R,2))+pow(2*math.pi*a/wavelength,2))*(math.cos(2*math.pi*R/wavelength)-j*math.sin(2*math.pi*R/wavelength))/(8*pow(math.pi,2)*wavelength*pow(R/wavelength,3)))

for m in range(0,subdivisions):
	for n in range(0,subdivisions):

		zm = Zs[m]
		ZN = Zs[n]
		Ir = sp.integrate.quad(real_eq2,ZN-(deltaZ/2),ZN+(deltaZ/2),args=(zm))
		Ii = sp.integrate.quad(imag_eq2,ZN-(deltaZ/2),ZN+(deltaZ/2),args=(zm))

		Z[m,n] = complex(Ir[0],Ii[0])*complex(0.0,integrand_C)

print Z

#Z = linalg.inv(Z)

I = linalg.solve(Z,V)


IFinal = []
EFinal_theta = []
Norm_Rad = []

for q in range(0,subdivisions):
	IFinal.append(abs(I[q,0]))
	summation_const = math.cos((k*antenna_length/2.0)+(k*deltaZ*Zs[q]))
	E[q,0] = j*eta*k*cmath.exp(-1.0*j*k*r)*I[q,0]*summation_const/(2.0*math.pi*r);
	print E[q,0]

for v in range(0,len(Angles)):
	tempsum = 0

	for q in range(0,subdivisions):
		#tempsum += E[q,0]*math.sin(Angles[v]) #for Etheta xz,yz cut planes
		tempsum += E[q,0]*math.sin(math.pi/2.0) #for Etheta xy cut plane


	tempsum = abs(tempsum)

	if tempsum == 0.0:
		EFinal_theta.append(0.0)
	else:
		EFinal_theta.append(10.0*math.log10(tempsum))

	Norm_Rad.append(abs(tempsum*tempsum)) 

Zin = (1/deltaZ)/IFinal[middle-1]

print "Input Impedance = " + str(Zin)

#plt.plot(Zs,IFinal)
#plt.ylabel('Current Amperes')
#plt.xlabel('Wire (distance m)')

plt.polar(Angles,Norm_Rad)
plt.title('Etheta Z = 0, xy-plane')
plt.show()





