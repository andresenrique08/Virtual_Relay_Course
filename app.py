# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 23:10:15 2020

@author: kraftmann
"""
import matplotlib.pyplot as plt
import pyComtrade
import scipy
from scipy.fftpack import fft, fftfreq, ifft
import scipy.signal as signal
import numpy as np
from tkinter import *
from tkinter import filedialog

# app
#funcion para cargar

comtradeObj = pyComtrade.ComtradeRecord()
def abrir(comtradeObj):
    ar=filedialog.askopenfilename(title="Abrir")#,filetypes=(("Archivos COMTRADE .cfg y .dat",".")))
   
    dat=filedialog.askopenfilename(title="Abrir")#,filetypes=(("Archivos COMTRADE.dat ",".dat")))
    
    comtradeObj.read(ar,dat)
    number_of_samples = comtradeObj['endsamp'][-1]
    sampling_freq=comtradeObj['samp'][-1]
    an=('Nombre de las señales analogas',comtradeObj.get_analog_ids())# print the ids of the analog channels.
    b=('Nombre de las señales digitales',comtradeObj.get_digital_ids())# print the ids of the analog channels.
    c=('Record has {} samples'.format(number_of_samples))
    d=('Sampling rate is {} samples/sec.'.format(comtradeObj['samp'][-1]))
    messagebox.showinfo("Información de los datos" ,an)
    messagebox.showinfo("Información de los datos" ,b)
    messagebox.showinfo("Información de los datos" ,c)
    messagebox.showinfo("Información de los datos" ,d)

def inicial(comtradeObj):
    print(len(comtradeObj['A'][0]['values']))
    time = comtradeObj.get_timestamps()
    print(len(time))
    # Reading voltaje and currents
    voltages = np.empty(([len(time),3]))
    currents = np.empty(([len(time),3]))

    for i in np.arange(6):
        if i<3:
            voltages[:,i] = comtradeObj['A'][i]['values']
            
        else:
            currents[:,i-3] = comtradeObj['A'][i]['values']

    f, axarr = plt.subplots(1, 2, figsize =(16, 4))
    f.suptitle('Lectura del archivo COMTRADE', y=1, fontsize=16)

    axarr[0].plot(time, voltages[:,0], 'b-', label='Phase A')
    axarr[0].plot(time, voltages[:,1], 'r-', label='Phase B')
    axarr[0].plot(time, voltages[:,2], 'g-', label='Phase C')
    axarr[0].set_xlabel('Time [sec]')
    axarr[0].set_ylabel('Voltage [V]')
    axarr[0].grid()
    axarr[0].legend()

    axarr[1].plot(time, currents[:,0], 'b-', label='Phase A')
    axarr[1].plot(time, currents[:,1], 'r-', label='Phase B')
    axarr[1].plot(time, currents[:,2], 'g-', label='Phase C')
    axarr[1].set_xlabel('Time [sec]')
    axarr[1].set_ylabel('Current [A]')
    axarr[1].grid()
    axarr[1].legend()

def ayuda():
    messagebox.showinfo("Carga de Archivos","cargue primero el archivo .cfg y despues el .dat")

#fs_user_cycle=StringVar()

def guardara():
    e1=en1.get()
    e2=en2.get()
    if e1=="" and e2=="":
        messagebox.showinfo("Introduzca parametros", "Introduzca la frecuencia de muestreo del rele y la relación de transformación" )
    
    
    
def subsampling(time,data,fs_comtrade,fk,fs_user_cycle):
    # time is the vector of time
    # data is the vector with the signal
    # fs_comtrade is the sample rate from the comtrade file
    # fk is the frequency of the system
    # fs_user_cycle is the sample rate given by user
    N1 = fs_user_cycle
    fs_cycle = fs_comtrade/fk
    N=np.int(fs_cycle)
    N_tot = np.int(len(data)/fs_cycle)
    new_data = [0]
    new_time = [0]
    for i in np.arange(N_tot):
        xi=data[i*N:i*N+N]
        ti=time[i*N:i*N+N]
        new_data[i*N1:i*N1+N1] = signal.resample(xi, N1)
        new_time[i*N1:i*N1+N1] = np.linspace(ti[0], ti[-1], N1, endpoint=False)
        
    return (new_time,new_data)

    
raiz=Tk()
raiz.title("mi app")
#raiz.geometry("650x380")
mifr= Frame()
mifr.pack()
mifr.config()
mifr.config(width="650",height="380")
mifr.config(cursor="star")

bmenu=Menu(raiz)
raiz.config(menu=bmenu)

barchiv=Menu(bmenu,tearoff=0)
bhelp=Menu(bmenu,tearoff=0)
barchiv.add_command(label="Abrir archivo ", command=abrir)
barchiv.add_command(label="Guardar")
bmenu.add_cascade(label="Archivo",menu=barchiv)
#help

bhelp.add_command(label="Ayuda", command=ayuda)
bmenu.add_cascade(label="Paso a paso",menu=bhelp)


Button(raiz,text="Abrir archivos",command=abrir(comtradeObj)).pack()

bi=Button(mifr,text="Mostrar señales de entrada",command=abrir(comtradeObj))
bi.grid(row=4,column=1,sticky="w")

t1=Label(mifr,text="introduzca frecuencia de muestro del rele",fg="green")
t1.grid(row=0,column=0,sticky="w",pady="10")

en1=Entry(mifr)
en1.grid(row=0,column=1,sticky="w",pady="10")
b1=Button(mifr,text="cargar valores",command=guardara)
b1.grid(row=3,column=1,sticky="w")

t2=Label(mifr,text="introduzca Relación de transformadores de instrumentación ",fg="blue")
t2.grid(row=1,column=0,sticky="w")
en2=Entry(mifr)
en2.grid(row=1,column=1,sticky="w",pady="10")


raiz.mainloop()
