from scipy.fftpack import fft, fftfreq
import scipy.io.wavfile as wav
import numpy as np
import matplotlib.pylab as plt
#----------------------------------------------------------------
dat=wav.read("violin.wav")[1]# tomar los datos de un audio
sr=wav.read("violin.wav")[0] #velocidad de muestreo

n=len(dat)
print n
dt = 1. / ( sr) 
fft_x = fft(dat) / n # FFT Normalized
freq = fftfreq(n, dt) # Recuperamos las frecuencias
fo=abs(fft_x)
def filtro(f,t):
	for i in range(len(f)):
		if 1000>abs(f[i]) or abs(f[i]) >2000:
			t[i]=0.
	
		
	return t

plt.subplot(2,1,1)
plt.plot(freq,fo)
plt.subplot(2,1,2)
plt.plot(freq,filtro(freq,fo),'g')
plt.show()
	
