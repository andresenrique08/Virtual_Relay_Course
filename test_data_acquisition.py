# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 18:31:15 2020

@author: anavr
"""

import numpy as np
from scipy.signal import butter, lfilter, freqz
from matplotlib import pyplot as plt
import pyComtrade
comtradeObj = pyComtrade.ComtradeRecord()
comtradeObj.read('./Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.cfg', './Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.dat')

print('Los analogos',comtradeObj.get_analog_ids())  # print the ids of the analog channels.

N = comtradeObj['endsamp'][-1]

print('Record has {} samples'.format(N))
print('Sampling rate is {} samples/sec.'.format(comtradeObj['samp'][-1]))

# Reading Phase A:
Voltage_A = comtradeObj['A'][1]['values']
Current_A = comtradeObj['A'][4]['values']

# Reading time vector:
time = comtradeObj.get_timestamps()

# Defining methods for the filter
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=True)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order)
    y= lfilter(b,a,data)
    return y

# Filter requirements.
order = 2
fs = comtradeObj['samp'][-1]       # sample rate, Hz
cutoff = 60  # desired cutoff frequency of the filter, Hz

# Filter the data, and plot both the original and filtered signals.
VA_filter = butter_lowpass_filter(Voltage_A, cutoff, fs, order)
IA_filter = butter_lowpass_filter(Current_A, cutoff, fs, order)

plt.subplot(2, 1, 1)
plt.plot(time, Voltage_A, 'b-', label='data')
plt.plot(time, VA_filter, 'g-', linewidth=2, label='filtered data')
plt.xlabel('Time [sec]')
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(time, Current_A, 'b-', label='data')
plt.plot(time, IA_filter, 'g-', linewidth=2, label='filtered data')
plt.xlabel('Time [sec]')
plt.grid()
plt.legend()

plt.subplots_adjust(hspace=0.35)
plt.show()