#!/usr/bin/python

import numpy as np
import scipy as sp
from scipy.integrate import quad
from scipy import linalg
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

eta = 377.0
j = complex(0.0,1.0)
Io = 1.0
freq = 2.4e09
c = 3.0e8
lam = c/freq


k = math.pi*2.0/lam
length = lam/2.0
subdivisions = 100
segments = length/subdivisions
angleSteps = 240

V = np.zeros(shape=(angleSteps,subdivisions)) #solution vector
V = np.complex128(V)

PRad = []
PRad2 = []
PRadNew = []
PRadNew2 = []
field = []
field2 = []
PRadBent = []
Angles = []
Angles2 = []
alphaTheta = []

phaseLength = 0.0
summation2 = 0.0
summation = 0.0

#theta equations
def real_eq2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(-1.0*math.sin(theta)*math.cos(k*z)*cmath.exp(j*k*L))

def imag_eq2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(-1.0*math.sin(theta)*math.cos(k*z)*cmath.exp(j*k*L))

def real_eq2Z2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(-1.0*math.sin(theta)*math.cos(k*(z + length/4.0))*cmath.exp(j*k*L))

def imag_eq2Z2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(-1.0*math.sin(theta)*math.cos(k*(z + length/4.0))*cmath.exp(j*k*L))

def real_eq2X1(x,theta,z,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(math.cos(theta)*math.cos(phi)*math.cos(k*(x - length/4.0))*cmath.exp(j*k*L))

def imag_eq2X1(x,theta,z,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(math.cos(theta)*math.cos(phi)*math.cos(k*(x - length/4.0))*cmath.exp(j*k*L))

def real_eq2Y1(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(math.cos(theta)*math.sin(phi)*math.cos(k*(y + length/8.0))*cmath.exp(j*k*L))

def imag_eq2Y1(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(math.cos(theta)*math.sin(phi)*math.cos(k*(y + length/8.0))*cmath.exp(j*k*L))

def real_eq2Y2(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(math.cos(theta)*math.sin(phi)*math.cos(k*(length/2.0 - y))*cmath.exp(j*k*L))

def imag_eq2Y2(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(math.cos(theta)*math.sin(phi)*math.cos(k*(length/2.0 - y))*cmath.exp(j*k*L))
def bookeq(z,theta):
	return -2.0*math.sin(k*(length-z))*math.cos(k*z*math.cos(theta))

## APHI Equations #####
def PHIreal_eq2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(0.0*math.cos(k*z)*cmath.exp(j*k*L))

def PHIimag_eq2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(0.0*math.cos(k*z)*cmath.exp(j*k*L))

def PHIreal_eq2Z2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(0.0*math.cos(k*(z + length/4.0))*cmath.exp(j*k*L))

def PHIimag_eq2Z2(z,theta,x,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(0.0*math.cos(k*(z + length/4.0))*cmath.exp(j*k*L))

def PHIreal_eq2X1(x,theta,z,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(-1.0*math.sin(phi)*math.cos(k*(x - length/4.0))*cmath.exp(j*k*L))

def PHIimag_eq2X1(x,theta,z,y,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(-1.0*math.sin(phi)*math.cos(k*(x - length/4.0))*cmath.exp(j*k*L))

def PHIreal_eq2Y1(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(math.cos(phi)*math.cos(k*(y + length/8.0))*cmath.exp(j*k*L))

def PHIimag_eq2Y1(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(math.cos(phi)*math.cos(k*(y + length/8.0))*cmath.exp(j*k*L))

def PHIreal_eq2Y2(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.real(math.cos(phi)*math.cos(k*(length/2.0 - y))*cmath.exp(j*k*L))

def PHIimag_eq2Y2(y,theta,x,z,phi):
	L = x*math.sin(theta)*math.cos(phi) + y*math.sin(theta)*math.sin(phi) + z*math.cos(theta)
	#return sp.real(-1.0*math.sin(theta)*math.sin(k*(z + length/2.0))*cmath.exp(j*z*k*math.cos(theta)))
	return sp.imag(math.cos(phi)*math.cos(k*(length/2.0 - y))*cmath.exp(j*k*L))
##  END THETA EQUATIONS ###

#see page 84 in elliot
#for single dipole
for th in range(0,angleSteps):
	theta = math.pi*th/120.0
	Angles.append(theta)

	Ir = sp.integrate.quad(real_eq2,-length/2.0,length/2.0,args=(theta,0,0,0))
	Ii = sp.integrate.quad(imag_eq2,-length/2.0,length/2.0,args=(theta,0,0,0))
	Izr = sp.integrate.quad(bookeq,0,length,args=(theta))

	myval = complex(Ir[0],Ii[0])

	#print myval
	myrecord = (abs(myval))**2*math.sin(theta)
	print myrecord
	myrecord2 = (abs(Izr[0])**2)*math.sin(theta)

	field.append(abs(myval))

	#append the absolute value in order to get a positive magnitude despite the sin coefficient for theta
	PRadNew.append(abs(myrecord))
	PRadNew2.append(myrecord2)


#TODO: need to still do the effective aperture
Rr = sum(PRadNew)*k**2*eta/(16.0*math.pi)*(math.pi/240.0)*2
print "Radiation Resistance %f" % Rr
#print sum(PRadNew2)

#For x-z plane
plt.polar(Angles,field)
plt.title('xz-plane')
plt.show()

#For y-z plane
#plt.polar(Angles,field)
#plt.title('yz-plane')
#plt.show()

#For x-y plane, change theta in integration above to be 0.
#plt.polar(Angles,field)
#plt.title('xy-plane')
#plt.show()
#for ph in range(0,angleSteps):
#	phi = math.pi*ph/120.0
for th in range(0,angleSteps):
	theta = math.pi*th/120.0
	Angles2.append(theta)

	#x1r = sp.integrate.quad(real_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,0))
	#x1i = sp.integrate.quad(imag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,0))
	x1r = sp.integrate.quad(real_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,math.pi/2.0))
	x1i = sp.integrate.quad(imag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,math.pi/2.0))
	#x1r = sp.integrate.quad(real_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,phi))
	#x1i = sp.integrate.quad(imag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,phi))
	x1c = complex(x1r[0],x1i[0])

	#z1r = sp.integrate.quad(real_eq2,-length/4.0,length/8.0,args=(theta,0,0,0))
	#z1i = sp.integrate.quad(imag_eq2,-length/4.0,length/8.0,args=(theta,0,0,0))
	z1r = sp.integrate.quad(real_eq2,-length/4.0,length/8.0,args=(theta,0,0,math.pi/2.0))
	z1i = sp.integrate.quad(imag_eq2,-length/4.0,length/8.0,args=(theta,0,0,math.pi/2.0))
	#z1r = sp.integrate.quad(real_eq2,-length/4.0,length/8.0,args=(theta,0,0,phi))
	#z1i = sp.integrate.quad(imag_eq2,-length/4.0,length/8.0,args=(theta,0,0,phi))
	z1c = complex(z1r[0],z1i[0])

	#y1r = sp.integrate.quad(real_eq2,0,length/8.0,args=(theta,0,length/8.0,0))
	#y1i = sp.integrate.quad(imag_eq2,0,length/8.0,args=(theta,0,length/8.0,0))
	y1r = sp.integrate.quad(real_eq2,0,length/8.0,args=(theta,0,length/8.0,math.pi/2.0))
	y1i = sp.integrate.quad(imag_eq2,0,length/8.0,args=(theta,0,length/8.0,math.pi/2.0))
	#y1r = sp.integrate.quad(real_eq2,0,length/8.0,args=(theta,0,length/8.0,phi))
	#y1i = sp.integrate.quad(imag_eq2,0,length/8.0,args=(theta,0,length/8.0,phi))
	y1c = complex(y1r[0],y1i[0])

	#z2r = sp.integrate.quad(real_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,0))
	#z2i = sp.integrate.quad(imag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,0))
	z2r = sp.integrate.quad(real_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,math.pi/2.0))
	z2i = sp.integrate.quad(imag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,math.pi/2.0))
	#z2r = sp.integrate.quad(real_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,phi))
	#z2i = sp.integrate.quad(imag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,phi))
	z2c = complex(z2r[0],z2i[0])

	#y2r = sp.integrate.quad(real_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,0))
	#y2i = sp.integrate.quad(imag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,0))
	y2r = sp.integrate.quad(real_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,math.pi/2.0))
	y2i = sp.integrate.quad(imag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,math.pi/2.0))
	#y2r = sp.integrate.quad(real_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,phi))
	#y2i = sp.integrate.quad(imag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,phi))
	y2c = complex(y2r[0],y2i[0])

	#PHI EQUATIONS
	#PHIx1r = sp.integrate.quad(PHIreal_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,0))
	#PHIx1i = sp.integrate.quad(PHIimag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,0))
	PHIx1r = sp.integrate.quad(PHIreal_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,math.pi/2.0))
	PHIx1i = sp.integrate.quad(PHIimag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,math.pi/2.0))
	#PHIx1r = sp.integrate.quad(PHIreal_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,phi))
	#PHIx1i = sp.integrate.quad(PHIimag_eq2X1,-length/4.0,0,args=(theta,-length/4.0,0,phi))
	PHIx1c = complex(PHIx1r[0],PHIx1i[0])

	#PHIz1r = sp.integrate.quad(PHIreal_eq2,-length/4.0,length/8.0,args=(theta,0,0,0))
	#PHIz1i = sp.integrate.quad(PHIimag_eq2,-length/4.0,length/8.0,args=(theta,0,0,0))
	PHIz1r = sp.integrate.quad(PHIreal_eq2,-length/4.0,length/8.0,args=(theta,0,0,math.pi/2.0))
	PHIz1i = sp.integrate.quad(PHIimag_eq2,-length/4.0,length/8.0,args=(theta,0,0,math.pi/2.0))
	#PHIz1r = sp.integrate.quad(PHIreal_eq2,-length/4.0,length/8.0,args=(theta,0,0,phi))
	#PHIz1i = sp.integrate.quad(PHIimag_eq2,-length/4.0,length/8.0,args=(theta,0,0,phi))
	PHIz1c = complex(PHIz1r[0],PHIz1i[0])

	#PHIy1r = sp.integrate.quad(PHIreal_eq2,0,length/8.0,args=(theta,0,length/8.0,0))
	#PHIy1i = sp.integrate.quad(PHIimag_eq2,0,length/8.0,args=(theta,0,length/8.0,0))
	PHIy1r = sp.integrate.quad(PHIreal_eq2,0,length/8.0,args=(theta,0,length/8.0,math.pi/2.0))
	PHIy1i = sp.integrate.quad(PHIimag_eq2,0,length/8.0,args=(theta,0,length/8.0,math.pi/2.0))
	#PHIy1r = sp.integrate.quad(PHIreal_eq2,0,length/8.0,args=(theta,0,length/8.0,phi))
	#PHIy1i = sp.integrate.quad(PHIimag_eq2,0,length/8.0,args=(theta,0,length/8.0,phi))
	PHIy1c = complex(PHIy1r[0],PHIy1i[0])

	#PHIz2r = sp.integrate.quad(PHIreal_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,0))
	#PHIz2i = sp.integrate.quad(PHIimag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,0))
	PHIz2r = sp.integrate.quad(PHIreal_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,math.pi/2.0))
	PHIz2i = sp.integrate.quad(PHIimag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,math.pi/2.0))
	#PHIz2r = sp.integrate.quad(PHIreal_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,phi))
	#PHIz2i = sp.integrate.quad(PHIimag_eq2Z2,length/8.0,length/4.0,args=(theta,0,length/8.0,phi))
	PHIz2c = complex(PHIz2r[0],PHIz2i[0])

	#PHIy2r = sp.integrate.quad(PHIreal_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,0))
	#PHIy2i = sp.integrate.quad(PHIimag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,0))
	PHIy2r = sp.integrate.quad(PHIreal_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,math.pi/2.0))
	PHIy2i = sp.integrate.quad(PHIimag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,math.pi/2.0))
	#PHIy2r = sp.integrate.quad(PHIreal_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,phi))
	#PHIy2i = sp.integrate.quad(PHIimag_eq2Y2,length/8.0,0,args=(theta,0,length/4.0,phi))
	PHIy2c = complex(PHIy2r[0],PHIy2i[0])

	comboval = x1c + z1c + y1c + z2c + y2c + PHIx1c + PHIz1c + PHIy1c + PHIz2c + PHIy2c
	#print myval
	myrecord = (abs(comboval))**2*math.sin(theta)

	field2.append(abs(comboval))

	PRadBent.append(abs(myrecord))


print sum(PRadBent)*k**2*eta/(16.0*math.pi)*(math.pi/240.0)
#print sum(PRadNew2)

#xz-plane
#plt.polar(Angles2,PRadBent)
#plt.title('xz-plane')
#plt.show()

#yz-plane
plt.polar(Angles2,PRadBent)
plt.title('yz-plane')
plt.show()

#xy-plane enable the phi'd equations for this...
#plt.polar(Angles2,PRadBent)
#plt.title('yz-plane')
#plt.show()

#yz plane change phi to 90






