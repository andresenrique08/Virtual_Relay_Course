# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:30:50 2020

@author: anavr
"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

order = 10
normal_cutoff = 15
b, a = signal.butter(order, normal_cutoff, 'low', analog=True)
w, h = signal.freqs(b, a)
#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

plt.subplot(2,1,1)
plt.title('Butterworth filter frequency response')
plt.semilogx(w, 20 * np.log10(abs(h)))
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(normal_cutoff, color='green') # cutoff frequency
plt.tight_layout()



t = np.linspace(0, 1, 1000, False)  # 1 second
sig = np.sin(2*np.pi*10*t) + np.sin(2*np.pi*20*t)

sos = signal.butter(10, 15, 'hp', fs=1000, output='sos')
filtered = signal.sosfilt(sos, sig)
#filtered = signal.lfilter(b,a, sig)

plt.subplot(2,1,2)
plt.plot(t, sig, label='10 Hz and 20 Hz sinusoids')
plt.plot(t, filtered, label='After 15 Hz high-pass filter')
plt.axis([0, 1, -2, 2])
plt.xlabel('Time [seconds]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
