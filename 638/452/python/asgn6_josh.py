import numpy as np
from scipy import linalg
from scipy.integrate import quad
from scipy import real
from scipy import imag
from matplotlib import pyplot
from cmath import cos
from cmath import sin
#CONSTANTS
L = 0.47		#wavelength = 1 meter
N = 11
deltaz = L/N
#print ("deltaz=")
#print (deltaz)
invdeltaz = N/L
a = 0.005
eo = np.pi * 4e-7
C = (1/(1j))*(1/(2 * np.pi * 300e6 * eo))*(1/(4*np.pi))
#print ("C=")
#print (C)
k = 2 * np.pi / 1

#CREATE VOLTAGE VECTOR
Vn = np.zeros( (N,1), dtype = complex )
feed = ((N+1)/2) - 1
Vn[feed] = invdeltaz
#print ("Vn=")
print (Vn)

#CREATE Zn and Zm VECTORS FOR IMPEDANCE CALCULATIONS
Zn = np.zeros((N,1))
for i in range(N):
	Zn[i][0]=-1*(((N+1)/2)-1)*deltaz + i*deltaz
Zm=Zn
#print (Zm)
#print ("Zn=")
#print (Zn)

#FUNCTION TO CALCULATE COMPLEX INTEGRAL (TAKEN/MODIFIED FROM INTERNET)
def complex_quadrature(func, a, b, **kwargs):
	def real_func(x):
		return real(func(x))
	def imag_func(x):
		return imag(func(x))
	real_integral = quad(real_func, a, b, **kwargs)
	imag_integral = quad(imag_func, a, b, **kwargs)
	return (real_integral[0] + 1j*imag_integral[0])

#CREATE IMPEDANCE MATRIX
Z=np.zeros( (N,N), dtype=complex) #CREATE Z
#print (Z)
for n in range(N):
	for m in range(N):
		def Zmn(Zprime): #Zmn DEFINED IN TERMS OF Zprime
#		    R = (np.sqrt(a**2+(Zm[m]-Zprime)**2))
		    return (C*((np.cos(1j*k*(np.sqrt(a**2+(Zm[m]-Zprime)**2)))-1j*np.sin(1j*k*(np.sqrt(a**2+(Zm[m]-Zprime)**2))))/(np.sqrt(a**2+(Zm[m]-Zprime)**2))**5)*((1+1j*k*(np.sqrt(a**2+(Zm[m]-Zprime)**2)))*(2*(np.sqrt(a**2+(Zm[m]-Zprime)**2))**2-3*a**2)+(k*a*(np.sqrt(a**2+(Zm[m]-Zprime)**2)))**2))
#			return (C*(cos(1j*k*(np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 )))-1j*sin(1j*k*(np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 ))*(1/((np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 )) ** 5))*((1+1j*k*(np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 )))*(2*((np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 ))**2) - 3*(a**2)) + (k*a*(np.sqrt( a**2 + (Zm[m]-Zprime) ** 2 )))**2))
		result = complex_quadrature(Zmn, (Zn[n]-deltaz/2), (Zn[n]+deltaz/2))
		Z[n][m] = result	#PUT RESULT OF INTEGRAL INTO MATRIX
#print ("Z=")
print (Z)
ZINVERSE = linalg.inv(Z) #CALCULATE Z INVERSE
#print ("ZINVERSE=")
#print (ZINVERSE)

#CREATE CURRENT VECTOR
I=np.zeros( (N,1), dtype=complex)
I = ZINVERSE.dot(Vn)  #CURRENT IS RESULT OF MULTIPLICATION OF ZINVERSE AND Vn
#print ("I=")
#print (I)
I_magnitude=np.zeros( (N,1))	#FIND MAGNITUDE OF CURRENT FOR PLOTTING
for i in range(N):    #MAGNITUDE OF A+Bj = SQRT(A^2+B^2)
	I_magnitude[i][0]=np.sqrt( real(I[i])**2 + imag(I[i])**2)
#print ("I_magnitude=")
print (I_magnitude)

#INPUT IMPEDANCE AT FEED POINT:
Zin = Vn[feed]/I[feed]
#print ("Zin=")
print (Zin)

#PLOT MAGNITUDE OF CURRENT
for i in range(N):
	pyplot.plot(Zn[i],I_magnitude[i], 'bs')
pyplot.ylabel('I_magnitude')
pyplot.xlabel('Zn')
pyplot.show()
