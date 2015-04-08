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

freq = 2.45e09

lengths = (kiloHigh - kiloLow)/subdivisions


airVP = 3.0e8 #m
rg60VP = 1.988e8 # m/s

#WR340
widthA = 0.086 #for WR340 10cm
fmn = (airVP/2.0) * (1/0.086)
print fmn/freq
wr340VP = airVP * math.sqrt(1 - pow((fmn/freq),2))

#rogers
microstripVP = airVP/math.sqrt(2.2)


air = []
rg60 = []
wav340 = []
microstrip = []
xaxis = []

for q in range(0,subdivisions):
	air.append( (kiloLow + q*lengths)/airVP)
	rg60.append( (kiloLow + q*lengths)/rg60VP)
	wav340.append( (kiloLow + q*lengths)/wr340VP)
	microstrip.append((kiloLow + q*lengths)/microstripVP)
	xaxis.append((kiloLow + q*lengths))



line1 = plt.plot(xaxis,air, 'b', label = 'Air')
line2 = plt.plot(xaxis,rg60, 'r', label = 'rg60')
line4 = plt.plot(xaxis,wav340, 'y', label = 'WR340')
line6 = plt.plot(xaxis,microstrip, 'c', label = 'microstrip')

#plt.legend([line1, line2, line3, line4, line5, line6], ['Air', 'RG60', 'CircWav', 'WR340', 'Fiber', 'Microstrip'])
#red_patch = mpatches.Patch(color='red', label='The red data')
#plt.legend(handles=[red_patch])

#plt.yscale('log')

legend = plt.legend(loc='upper center', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

plt.ylabel('Seconds')
plt.xlabel('km')

plt.show()








