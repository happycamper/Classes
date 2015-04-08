#!/usr/bin/python

import numpy as np
import scipy as sp
from scipy.integrate import quad
from scipy import linalg
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



freq = 2.45e09
c = 3e8
lam = c/freq
k = 2*math.pi/lam
l = 0.1 * lam
radCoeff = (1-math.cos(k*l/2.0))
j = complex(0.0,1.0)
eta = 377.0
IA = 1.0



subdivisions = 10000
kiloLow = .1*lam
kiloHigh = 100*lam


lengths = (kiloHigh - kiloLow)/subdivisions

wFar = []
wTot = []
xaxis = []


for q in range(0,subdivisions):
	dist = (kiloLow + q*lengths)
	wFar.append( (eta*IA**2*radCoeff**2) / (8*math.pi**2*dist**2))
	magImag = abs(complex((eta*IA**2*radCoeff**2) / (8*math.pi**2*dist**2),(eta/k)*(IA*.1*lam/(4*math.pi))**2*(1/pow(dist,5))))
	wTot.append(magImag)
	xaxis.append(dist/lam)



line1 = plt.plot(xaxis,wFar, 'b', label = 'WFarField')
line2 = plt.plot(xaxis,wTot, 'r', label = 'WTotal')


#plt.legend([line1, line2, line3, line4, line5, line6], ['Air', 'RG60', 'CircWav', 'WR340', 'Fiber', 'Microstrip'])
#red_patch = mpatches.Patch(color='red', label='The red data')
#plt.legend(handles=[red_patch])

#plt.yscale('log')
plt.xscale('log')

legend = plt.legend(loc='upper center', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

plt.ylabel('W/m^3')
plt.xlabel('wavelengths')

plt.show()








