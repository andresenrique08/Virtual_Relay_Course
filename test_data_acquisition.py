# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 18:31:15 2020

@author: anavr
"""

import numpy as np
from scipy import signal
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
def butter_lowpass(cutoff, fs, order):
    nyq = fs/2
    normal_cutoff = cutoff/nyq
    b,a = signal.butter(order, normal_cutoff, 'low', analog=True)
    return b,a

cutoff = 50
fs = np.ceil(1/cutoff * comtradeObj['samp'][-1])
# fs = comtradeObj['samp'][-1]/50
order = 2
b,a = butter_lowpass(cutoff, fs, order)

w, h = signal.freqs(b, a)
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# plt.subplot(2,1,1)
# plt.title('Butterworth filter frequency response')
# #plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
# plt.semilogx(0.5*fs*w/np.pi, 20 * np.log10(abs(h)))
# plt.xlabel('Frequency [radians / second]')
# plt.ylabel('Amplitude [dB]')
# plt.margins(0, 0.1)
# plt.grid(which='both', axis='both')
# plt.axvline(cutoff, color='green') # cutoff frequency
# plt.tight_layout()

zi = signal.lfilter_zi(b, a)
VA_Filtered = signal.lfilter(b,a,Voltage_A)
IA_Filtered = signal.lfilter(b,a,Current_A)


# print(VA_Filtered)
# print(len(Voltage_A))

plt.subplot(2, 1, 1)
plt.plot(time, Voltage_A, 'b-', label='data')
plt.plot(time, VA_Filtered, 'g-', linewidth=2, label='filtered data')
plt.xlabel('Time [sec]')
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(time, Current_A, 'b-', label='data')
plt.plot(time, IA_Filtered, 'g-', linewidth=2, label='filtered data')
plt.xlabel('Time [sec]')
plt.grid()
plt.legend()

# plt.subplots_adjust(hspace=0.35)
# plt.show()