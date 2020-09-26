# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:06:16 2020

@author: anavr
"""

import matplotlib.pyplot as plt
import pyComtrade
import scipy
from scipy.fftpack import fft, fftfreq, ifft
import scipy.signal as signal
import numpy as np
from tkinter import *
from tkinter import filedialog

comtradeObj = pyComtrade.ComtradeRecord()
comtradeObj.read('./Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.cfg', './Comtrade TEST 1-1/L1L2L3_75%/BINAIRE/Test_1-1_ABC60_75%_BIN.dat')


print('Los analogos',comtradeObj.get_analog_ids())  # print the ids of the analog channels.
print('Los digitales',comtradeObj.get_digital_ids())  # print the ids of the analog channels.

number_of_samples = comtradeObj['endsamp'][-1]
sampling_freq=comtradeObj['samp'][-1]

time = comtradeObj.get_timestamps() 

print('Record has {} samples'.format(number_of_samples))
print('uerSampling rate is {} samples/sec.'.format(comtradeObj['samp'][-1]))

 

cutoff_freq = 2400 # 60 hertz
sampling_duration = number_of_samples / sampling_freq # sampling duration in [s]
W= 2 * cutoff_freq / sampling_freq # Normalize the frequency


# Reading channel 4:
Va = comtradeObj['A'][0]['values']

fig, ax = plt.subplots()
ax.plot(time,Va,'r')
ax.set(xlabel='time (s)', ylabel=' Magnitud voltage ',
       title='phase A B C (V)')
ax.grid()

#plt.show()

Vb = comtradeObj['A'][1]['values']
plt.plot(time,Vb,'b')
#plt.show()

Vc = comtradeObj['A'][2]['values']
plt.plot(time,Vc,'g')
plt.show()

fig, ax = plt.subplots()
ax.set(xlabel='time (s)', ylabel=' Magnitud Corriente ',
       title='phase A B C (A)')
ax.grid()
Ia = comtradeObj['A'][3]['values']
plt.plot(time,Ia,'r')
#plt.show()

Ib = comtradeObj['A'][4]['values']
plt.plot(time,Ib,'b')
#plt.show()

Ic = comtradeObj['A'][5]['values']
plt.plot(time,Ic,'g')
plt.show()



    

    
fmr=32 #freceuncia de muestreo del rele, user define


# Filtros 
#butterlow


b,a= scipy.signal.butter(2, W, 'low', analog = True)
filtered_signal = scipy.signal.lfilter(b, a, Va )

V=np.array([[Va], [Vb],[Vc]])

If=np.array([[Ia],[Ib],[Ic]])


f,c,n=np.shape(V)

vt=np.transpose(Va)

Vf=np.empty([f, np.size(Va)])
Iff=np.empty([f, np.size(Ia)])
vn=range(len(V))
col=['r','g','b']
ti=['Señal filtrada fase A','Señal filtrada fase B','Señal filtrada fase C']

# Voltajes fistrados
for i in vn:
    
    Vf[i,:]=scipy.signal.lfilter(b, a, V[i,:])
    
    plt.plot(time, Vf[i,:], col[i])
    plt.ylabel('magnitud Voltaje')
    plt.xlabel('Time [s]')
    plt.grid
    plt.title(ti[i])
    #ax.set(xlabel='time (s)', ylabel=' Magnitud voltage ',titel=ti[i])
    plt.show()
    
# Corrientes filtradas

    
    Iff[i,:]=scipy.signal.lfilter(b, a, If[i,:])
    
    plt.plot(time, Iff[i,:], col[i])
    plt.ylabel('magnitud corriente')
    plt.xlabel('Time [s]')
    plt.grid
    plt.title(ti[i])
    plt.show()
        




dt=1/sampling_freq
fft_filt=fft(filtered_signal)
ffilt=abs(fft_filt)
freq1=fftfreq(n,dt)
plt.plot(freq1,ffilt,'g-',)
plt.show()



