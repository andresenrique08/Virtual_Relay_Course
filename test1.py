# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:06:16 2020

@author: anavr
"""

import matplotlib.pyplot as plt
import pyComtrade
comtradeObj = pyComtrade.ComtradeRecord()
comtradeObj.read('./Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.cfg', './Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.dat')

print(comtradeObj.get_analog_ids())  # print the ids of the analog channels.

N = comtradeObj['endsamp'][-1]

print('Record has {} samples'.format(N))
print('Sampling rate is {} samples/sec.'.format(comtradeObj['samp'][-1]))

# Reading channel 4:
AnalogChannelData = comtradeObj['A'][22]['values']

DigitalChannelData = comtradeObj['D'][25]['values']

# Reading time vector:
time = comtradeObj.get_timestamps()

# Ploting with matplotlib
# pylab.plot(time,channelData)
f, axarr = plt.subplots(2, sharex=True)

axarr[0].plot(time, AnalogChannelData)
axarr[0].set_title('pyComtrade Demo')
axarr[0].grid()
axarr[1].plot(time, DigitalChannelData)
axarr[1].set_ylim(top=1.05)  # bottom unchanged
axarr[1].grid()
axarr[1].set_xlabel('Time [s]')
plt.show()