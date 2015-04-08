#!/usr/bin/python

import numpy as np
import scipy as sp
from scipy.integrate import quad
from scipy import linalg
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

subdivisions = 10000
kiloLow = 1*1000 #m
kiloHigh = 10000*1000

lengths = (kiloHigh - kiloLow)/subdivisions


attAir = 0.05/1000.0 #db/m
attRG60 = 0.81 #dB/m
circWavatt = 0.001 #dB/m
wr_340Att = 0.0018   #TE waveguide mode found on pg 115 of pozar
attFiber = 0.35/1000.0 #db/m


#microstrip stuff
mWidth = 2.68/1000.0 #m
mResist = 1.67e-06 #m
mZo = 50.0
attMicro = (8.69*mResist/(mWidth*mZo)) #add neper conversion

air = []
rg60 = []
circwav = []
wav340 = []
fiber = []
microstrip = []
xaxis = []

for q in range(0,subdivisions):
	air.append( (kiloLow + q*lengths)*attAir)
	rg60.append( (kiloLow + q*lengths)*attRG60)
	circwav.append( (kiloLow + q*lengths)*circWavatt)
	wav340.append( (kiloLow + q*lengths)*wr_340Att)
	fiber.append( (kiloLow + q*lengths)*attFiber)
	microstrip.append((kiloLow + q*lengths)*attMicro)
	xaxis.append((kiloLow + q*lengths))

for q in range(0,subdivisions):
	air[q] -= air[0]
	rg60[q] -= rg60[0]
	circwav[q] -= circwav[0]
	wav340[q] -= wav340[0]
	fiber[q] -= fiber[0]
	microstrip[q] -= microstrip[0]

#fig = plt.figure()
#ax = fig.add_subplot(2,1,1)

line1 = plt.plot(xaxis,air, 'b', label = 'Air')
line2 = plt.plot(xaxis,rg60, 'r', label = 'rg60')
line3 = plt.plot(xaxis,circwav, 'g', label = 'Circwav')
line4 = plt.plot(xaxis,wav340, 'y', label = 'WR340')
line5 = plt.plot(xaxis,fiber, 'k', label = 'Fiber')
line6 = plt.plot(xaxis,microstrip, 'c', label = 'microstrip')

#plt.legend([line1, line2, line3, line4, line5, line6], ['Air', 'RG60', 'CircWav', 'WR340', 'Fiber', 'Microstrip'])
#red_patch = mpatches.Patch(color='red', label='The red data')
#plt.legend(handles=[red_patch])

plt.yscale('log')
#plt.xscale('log')

legend = plt.legend(loc='upper center', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

plt.ylabel('Att (dB/m)')
plt.xlabel('km')

plt.show()








