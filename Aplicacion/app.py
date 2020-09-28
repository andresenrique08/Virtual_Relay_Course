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
from inspect import *


# app
#funcion para cargar

global comtradeObj 
comtradeObj = pyComtrade.ComtradeRecord()
         
def abrir():
    
    global fs_comtrade,fk,time,fs,voltages,currents,N
    
    ar=filedialog.askopenfilename(title="Abrir")#,filetypes=(("Archivos COMTRADE .cfg y .dat",".")))  
    dat=filedialog.askopenfilename(title="Abrir")#,filetypes=(("Archivos COMTRADE.dat ",".dat")))
    comtradeObj.read(ar,dat)
    N = comtradeObj['endsamp'][-1]
    #sampling_freq=comtradeObj['samp'][-1]
    fs_comtrade=comtradeObj['samp'][-1]
    an=('Nombre de las señales analogas',comtradeObj.get_analog_ids())# print the ids of the analog channels.
    b=('Nombre de las señales digitales',comtradeObj.get_digital_ids())# print the ids of the analog channels.
    c=('Record has {} samples'.format(fs_comtrade))
    d=('Sampling rate is {} samples/sec.'.format(comtradeObj['samp'][-1]))
    fk = comtradeObj['line_freq']
    time = comtradeObj.get_timestamps()
    voltages = np.empty(([len(time),3]))
    currents = np.empty(([len(time),3]))
    # Reading voltaje and currents
    for i in np.arange(6):
        if i<3:
            voltages[:,i] = comtradeObj['A'][i]['values']
            
        else:
            currents[:,i-3] = comtradeObj['A'][i]['values']
    #fs=
    #messagebox.showinfo("Información de los datos" ,an)
    #messagebox.showinfo("Información de los datos" ,b)
    #messagebox.showinfo("Información de los datos" ,c)
    #messagebox.showinfo("Información de los datos" ,d)
    


def inicial():
    
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
    plt.show()
    
    

def ayuda():
    messagebox.showinfo("Carga de Archivos","cargue primero el archivo .cfg y despues el .dat")

#fs_user_cycle=StringVar()

def guardara():
    global e1,e2,e3
    
    e1=int(en1.get())
    e2=int(en2.get())
    e3=int(en3.get())
    if e2=="" and e3=="":
        messagebox.showinfo("Introduzca parametros", "Introduzca la relación de transformación" )
    elif e2=="" :
        messagebox.showinfo("Introduzca parametros", "Introduzca la relación de transformación del primario" )
    elif e3=="":
        messagebox.showinfo("Introduzca parametros", "Introduzca la relación de transformación del secundario" )
    
    
    
    
    
    
    
    
def subsampling(data):
    # time is the vector of time
    # data is the vector with the signal
    # fs_comtrade is the sample rate from the comtrade file
    # fk is the frequency of the system
    # fs_user_cycle is the sample rate given by user
    global N1,fs_cycle
    
   
    # time is the vector of time
    # data is the vector with the signal
    # fs_comtrade is the sample rate from the comtrade file
    # fk is the frequency of the system
    
    N1 = e1# fs_user_cycle=e1 is the sample rate given by user
    fs_cycle = fs_comtrade/fk
    N=np.int(fs_cycle)
    N_tot = np.int(len(data)/fs_cycle) # qe es esta monda 
    #N_tot = np.int((N*fk/fs_comtrade)*e1)
    new_data = [0]
    new_time = [0]
    
    for i in np.arange(N_tot):
        xi=data[i*N:i*N+N]
        ti=time[i*N:i*N+N]
        new_data[i*N1:i*N1+N1] = signal.resample(xi, N1)
        new_time[i*N1:i*N1+N1] = np.linspace(ti[0], ti[-1], N1, endpoint=False)
        
    return (new_time,new_data)
def plot():
    
    global V_sub,C_sub,time_sub
    N_tot = np.int(N*fk/fs_comtrade)*e1
    V_sub = np.empty(([N_tot,3]))
    C_sub = np.empty(([N_tot,3]))
    time_sub = np.empty(([N_tot,6]))
    
    for i in np.arange(6):
        if i<3:
            time_sub[:,i], V_sub[:,i] = subsampling(voltages[:,i])
        else:
            time_sub[:,i], C_sub[:,i-3] = subsampling(currents[:,i-3])
    
    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('Submuestreo de las señales', y=0.92, fontsize=16)
    
    # Plot Voltages
    axarr[0,0].plot(time, voltages[:,0], 'b-', label='Phase A')
    axarr[0,0].plot( time_sub[:,0], V_sub[:,0], 'co-', label='Phase A resampled')
    axarr[1,0].plot(time, voltages[:,1], 'r-', label='Phase B')
    axarr[1,0].plot( time_sub[:,1], V_sub[:,1], 'mo-', label='Phase B resampled')
    axarr[2,0].plot(time, voltages[:,1], 'g-', label='Phase C')
    axarr[2,0].plot( time_sub[:,2], V_sub[:,2], 'yo-', label='Phase C resampled')
    
    # Plot Currents
    axarr[0,1].plot(time, currents[:,0], 'b-', label='Phase A')
    axarr[0,1].plot( time_sub[:,3], C_sub[:,0], 'co-', label='Phase A resampled')
    axarr[1,1].plot(time, currents[:,1], 'r-', label='Phase B')
    axarr[1,1].plot( time_sub[:,4], C_sub[:,1], 'mo-', label='Phase B resampled')
    axarr[2,1].plot(time, currents[:,2], 'g-', label='Phase C')
    axarr[2,1].plot( time_sub[:,5], C_sub[:,2], 'yo-', label='Phase C resampled')
    
    for i in np.arange(3):
        axarr[i,0].set_xlabel('Time [sec]')
        axarr[i,0].set_ylabel('Voltage [V]')
        axarr[i,0].grid()
        axarr[i,0].legend()
        axarr[i,1].set_xlabel('Time [sec]')
        axarr[i,1].set_ylabel('Current [A]')
        axarr[i,1].grid()
        axarr[i,1].legend()
    plt.show()
    
def quantizer(data, quantizing_bits):
    # data is the vector with the signal
    # quantizing_bits is the number of bits for the converter
    # Quantizer - S&H and ADC
    quantizing_levels   = 2 ** quantizing_bits
    quantizing_step     = (np.max(data)-np.min(data)) / quantizing_levels
    quantizing_signal   = np.round (data / quantizing_step) * quantizing_step;
    
    return quantizing_signal

def plot3():
    quantizing_bits_V = 4 # Valor típico: 12 (Voltaje)
    quantizing_bits_C = 8 # Valor típico: 16 (Corriente)
    N_tot = np.int(N*fk/fs_comtrade)*e1
    global dig_V_sub,dig_C_sub,fs_user_cycle
    fs_user_cycle=e1
    dig_V_sub = np.empty(([N_tot,3]))
    dig_C_sub = np.empty(([N_tot,3]))
    for i in np.arange(6):
        if i<3:
            dig_V_sub[:,i] = quantizer(V_sub[:,i], quantizing_bits_V)
        else:
            dig_C_sub[:,i-3] = quantizer(C_sub[:,i-3], quantizing_bits_C)
    
    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('Digitalización de la señal', y=0.92, fontsize=16)
    
    # Plot Voltages
    axarr[0,0].plot(time_sub[:,0][0:fs_user_cycle], V_sub[:,0][0:fs_user_cycle], 'b-', label='Phase A')
    axarr[0,0].plot( time_sub[:,0][0:fs_user_cycle], dig_V_sub[:,0][0:fs_user_cycle], 'c-', label='Phase A digital')
    axarr[1,0].plot(time_sub[:,1][0:fs_user_cycle], V_sub[:,1][0:fs_user_cycle], 'r-', label='Phase B')
    axarr[1,0].plot( time_sub[:,1][0:fs_user_cycle], dig_V_sub[:,1][0:fs_user_cycle], 'm-', label='Phase B digital')
    axarr[2,0].plot(time_sub[:,2][0:fs_user_cycle], V_sub[:,2][0:fs_user_cycle], 'g-', label='Phase C')
    axarr[2,0].plot( time_sub[:,2][0:fs_user_cycle], dig_V_sub[:,2][0:fs_user_cycle], 'y-', label='Phase C digital')
    
    # Plot Currents
    axarr[0,1].plot(time_sub[:,3][0:fs_user_cycle], C_sub[:,0][0:fs_user_cycle], 'b-', label='Phase A')
    axarr[0,1].plot( time_sub[:,3][0:fs_user_cycle], dig_C_sub[:,0][0:fs_user_cycle], 'c-', label='Phase A digital')
    axarr[1,1].plot(time_sub[:,4][0:fs_user_cycle], C_sub[:,1][0:fs_user_cycle], 'r-', label='Phase B')
    axarr[1,1].plot( time_sub[:,4][0:fs_user_cycle], dig_C_sub[:,1][0:fs_user_cycle], 'm-', label='Phase B digital')
    axarr[2,1].plot(time_sub[:,5][0:fs_user_cycle], C_sub[:,2][0:fs_user_cycle], 'g-', label='Phase C')
    axarr[2,1].plot( time_sub[:,5][0:fs_user_cycle], dig_C_sub[:,2][0:fs_user_cycle], 'y-', label='Phase C digital')
    
    for i in np.arange(3):
        axarr[i,0].set_xlabel('Time [sec]')
        axarr[i,0].set_ylabel('Voltage [V]')
        axarr[i,0].grid()
        axarr[i,0].legend()
        axarr[i,1].set_xlabel('Time [sec]')
        axarr[i,1].set_ylabel('Current [A]')
        axarr[i,1].grid()
        axarr[i,1].legend()
    plt.show()



def walsh(data):
    # time is the vector of time
    # data is the vector with the signal
    # fk is the frequency of the system
    # fs_user_cycle is the sample rate given by user
    
    Nn=np.int(fs_user_cycle)
    N_tot =(len(data)-Nn)
    Xc = [0]*N_tot
    Xs = [0]*N_tot
    t = [0]*N_tot
    
    for i in np.arange(N_tot):
        xi=data[i:i+Nn]
        t[i]=time[i]
        Xc_sum = 0
        Xs_sum = 0
        for k in np.arange(Nn):
            Xc_temp=xi[k]*np.cos(2*np.pi*k/(Nn))
            Xc_sum=Xc_sum+Xc_temp
            Xs_temp=xi[k]*np.sin(2*np.pi*k/(Nn))
            Xs_sum=Xs_sum+Xs_temp
            
        Xc[i]= 2/(Nn)*Xc_sum
        Xs[i]= 2/(Nn)* Xs_sum
        
    return t, Xc, Xs

def plot4():
    #Nn=np.int(fs_user_cycle)
    N_tot=np.int(len(dig_V_sub[:,1]))
    
    global Y_C,X_C,Y_V,X_V,Xs_C,Xs_V,Xc_V,N_tot_walsh,Xc_C,t 
    
    N_tot_walsh = np.int(N_tot-fs_user_cycle)
    Xc_V = np.empty(([N_tot_walsh,3]))
    Xs_V = np.empty(([N_tot_walsh,3]))
    Xc_C = np.empty(([N_tot_walsh,3]))
    Xs_C = np.empty(([N_tot_walsh,3]))
    X_V = np.empty(([N_tot_walsh,3]))
    Y_V = np.empty(([N_tot_walsh,3]))
    X_C = np.empty(([N_tot_walsh,3]))
    Y_C = np.empty(([N_tot_walsh,3]))
    t = np.empty(([N_tot_walsh,6]))
    
    for i in np.arange(6):
        if i<3:
            t[:,i], Xc_V[:,i], Xs_V[:,i] = walsh( dig_V_sub[:,i])
            X_V[:,i] = np.sqrt(np.power(Xc_V[:,i],2)+np.power(Xs_V[:,i],2))
            Y_V[:,i] = np.arctan(Xs_V[:,i]/Xc_V[:,i])*180/np.pi
        else:
            t[:,i], Xc_C[:,i-3], Xs_C[:,i-3] = walsh(dig_C_sub[:,i-3], )
            X_C[:,i-3] = np.sqrt(np.power(Xc_C[:,i-3],2)+np.power(Xs_C[:,i-3],2))
            Y_C[:,i-3] = np.arctan(Xs_C[:,i-3]/Xc_C[:,i-3])*180/np.pi
    
    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('FFT En Magnitud', y=0.92, fontsize=16)
    
    # Plot Voltages
    axarr[0,0].plot(time_sub[:,0], V_sub[:,0], 'b-', label='Phase A')
    axarr[0,0].plot( t[:,0], X_V[:,0], 'c-', label='Phase A FFT(mag)')
    axarr[1,0].plot(time_sub[:,1], V_sub[:,1], 'r-', label='Phase B')
    axarr[1,0].plot( t[:,1], X_V[:,1], 'm-', label='Phase B FFT(mag)')
    axarr[2,0].plot(time_sub[:,2], V_sub[:,2], label='Phase C')
    axarr[2,0].plot( t[:,2], X_V[:,2], 'y-', label='Phase C FFT(mag)')
    
    # Plot Currents
    axarr[0,1].plot(time_sub[:,3], C_sub[:,0], 'b-', label='Phase A')
    axarr[0,1].plot( t[:,3], X_C[:,0], 'c-', label='Phase A FFT(mag)')
    axarr[1,1].plot(time_sub[:,4], C_sub[:,1], 'r-', label='Phase B')
    axarr[1,1].plot( t[:,4], X_C[:,1], 'm-', label='Phase B FFT(mag)')
    axarr[2,1].plot(time_sub[:,5], C_sub[:,2], 'g-', label='Phase C')
    axarr[2,1].plot( t[:,5], X_C[:,2], 'y-', label='Phase C FFT(mag)')
    
    for i in np.arange(3):
        axarr[i,0].set_xlabel('Time [sec]')
        axarr[i,0].set_ylabel('Voltage [V]')
        axarr[i,0].grid()
        axarr[i,0].legend()
        axarr[i,1].set_xlabel('Time [sec]')
        axarr[i,1].set_ylabel('Current [A]')
        axarr[i,1].grid()
        axarr[i,1].legend()    
    plt.show()

    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('FFT En Magnitud', y=0.92, fontsize=16)
    
    # Plot Voltages
    axarr[0,0].plot(time_sub[:,0], V_sub[:,0], 'b-', label='Phase A')
    axarr[0,0].plot( t[:,0], Y_V[:,0], 'c-', label='Phase A FFT[ang(rad)]')
    axarr[1,0].plot(time_sub[:,1], V_sub[:,1], 'r-', label='Phase B')
    axarr[1,0].plot( t[:,1], Y_V[:,1], 'm-', label='Phase B FFT[ang(rad)]')
    axarr[2,0].plot(time_sub[:,2], V_sub[:,2], label='Phase C')
    axarr[2,0].plot( t[:,2], Y_V[:,2], 'y-', label='Phase C FFT[ang(rad)]')
    
    # Plot Currents
    axarr[0,1].plot(time_sub[:,3], C_sub[:,0], 'b-', label='Phase A')
    axarr[0,1].plot( t[:,3], Y_C[:,0], 'c-', label='Phase A FFT[ang(rad)]')
    axarr[1,1].plot(time_sub[:,4], C_sub[:,1], 'r-', label='Phase B')
    axarr[1,1].plot( t[:,4], Y_C[:,1], 'm-', label='Phase B FFT[ang(rad)]')
    axarr[2,1].plot(time_sub[:,5], C_sub[:,2], 'g-', label='Phase C')
    axarr[2,1].plot( t[:,5], Y_C[:,2], 'y-', label='Phase C FFT[ang(rad)]')
    
    for i in np.arange(3):
        axarr[i,0].set_xlabel('Time [sec]')
        axarr[i,0].set_ylabel('Angle (°)')
        axarr[i,0].grid()
        axarr[i,0].legend()
        axarr[i,1].set_xlabel('Time [sec]')
        axarr[i,1].set_ylabel('Angle (°)')
        axarr[i,1].grid()
        axarr[i,1].legend()
    plt.show()

    
def valo():
     #Voltajes
    
    a=('Voltaje fase A:', "{:.2f}".format(np.double(X_V[-1:,0])), 'V', "{:.2f}".format(np.double(Y_V[-1:,0])),'°')
    b=('Voltaje fase B:', "{:.2f}".format(np.double(X_V[-1:,1])), 'V', "{:.2f}".format(np.double(Y_V[-1:,1])),'°')
    c=('Voltaje fase C:', "{:.2f}".format(np.double(X_V[-1:,2])), 'V', "{:.2f}".format(np.double(Y_V[-1:,2])),'°')
    
    # Corrientes
    d=('Corriente fase A:', "{:.2f}".format(np.double(X_C[-1:,0])), 'A', "{:.2f}".format(np.double(Y_C[-1:,0])),'°')
    e=('Corriente fase B:', "{:.2f}".format(np.double(X_C[-1:,1])), 'A', "{:.2f}".format(np.double(Y_C[-1:,1])),'°')
    f=('Corriente fase C:', "{:.2f}".format(np.double(X_C[-1:,2])), 'A', "{:.2f}".format(np.double(Y_C[-1:,2])),'°')
    messagebox.showinfo(a,b,c)
    messagebox.showinfo(d,e,f)

def fin():
    fig, ax = plt.subplots(1, 2, figsize =(16, 6))
    fig.suptitle('Diagrama fasorial de Voltaje y Corriente', y=0.95, fontsize=16)
    
    lim_axis_V = np.max([np.float(X_V[-1:,0]), np.float(X_V[-1:,1]),np.float(X_V[-1:,2])])
    lim_axis_C = np.max([np.float(X_C[-1:,0]), np.float(X_C[-1:,1]),np.float(X_C[-1:,2])])
    # Creando el punto de origen para los vectores
    x_pos = [0, 0,0] 
    y_pos = [0, 0,0]
    
    ax[0].quiver(x_pos, y_pos, Xc_V[-1,0], Xs_V[-1,0], angles='xy', scale_units = 'xy', scale=1, color=['b'], label='Fase A')  
    ax[0].quiver(x_pos, y_pos, Xc_V[-1,1], Xs_V[-1,1], angles='xy', scale_units = 'xy', scale=1, color=['r'], label='Fase B')  
    ax[0].quiver(x_pos, y_pos, Xc_V[-1,2], Xs_V[-1,2], angles='xy', scale_units = 'xy', scale=1, color=['g'], label='Fase C')  
    ax[0].axis([-1.2*lim_axis_V, 1.2*lim_axis_V, -1.2*lim_axis_V, 1.2*lim_axis_V]) 
    ax[0].set_title('Voltaje [V]')
    ax[0].legend() #<-- Se nombran las leyendas
    ax[0].grid(b=True, which='major') #<-- plot grid lines
    
    ax[1].quiver(x_pos, y_pos, Xc_C[-1,0], Xs_C[-1,0], angles='xy', scale_units = 'xy', scale=1, color=['b'], label='Fase A')  
    ax[1].quiver(x_pos, y_pos, Xc_C[-1,1], Xs_C[-1,1], angles='xy', scale_units = 'xy', scale=1, color=['r'], label='Fase B')  
    ax[1].quiver(x_pos, y_pos, Xc_C[-1,2], Xs_C[-1,2], angles='xy', scale_units = 'xy', scale=1, color=['g'], label='Fase C')  
    ax[1].axis([-1.2*lim_axis_C, 1.2*lim_axis_C, -1.2*lim_axis_C, 1.2*lim_axis_C]) 
    ax[1].set_title('Corriente [A]')
    ax[1].legend() #<-- Se nombran las leyendas
    ax[1].grid(b=True, which='major') #<-- plot grid lines
    plt.show()
#---------------------------------Interfaz Grafica con Tikinter--------------------------
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
barchiv.add_command(label="Abrir archivo ", command=abrir)#carga los archivos
barchiv.add_command(label="Guardar")
bmenu.add_cascade(label="Archivo",menu=barchiv)
#help

bhelp.add_command(label="Ayuda", command=ayuda)
bmenu.add_cascade(label="Paso a paso",menu=bhelp)


Button(raiz,text="Abrir archivos",command=abrir).pack()

bi=Button(mifr,text="Mostrar señales de entrada",command= inicial)#carga los archivos
bi.grid(row=4,column=1,sticky="w")

#img = ImageTk.PhotoImage(Image.open(path))
#panel = tk.Label(root, image = img)

t1=Label(mifr,text="introduzca frecuencia de muestro del rele",fg="green")
t1.grid(row=0,column=0,sticky="w",pady="10")

en1=IntVar()
d=OptionMenu(mifr, en1, "4","8","16","32")
d.grid(row=0,column=1,sticky="w",pady="10")
en1.set("4")
# Relación de transformadores de instrumentación le falta  dividir entre prim y secundario
t2=Label(mifr,text="introduzca Relación de transformadores de instrumentación ",fg="blue")
t2.grid(row=1,column=0,sticky="w")
en2=Entry(mifr)
en2.grid(row=1,column=1,sticky="w",pady="10")

en3=Entry(mifr)
en3.grid(row=1,column=2,sticky="w",pady="10")
# estapa de mostrar señal de entrada
#en1=Entry(mifr)
#en1.grid(row=0,column=1,sticky="w",pady="10")
b1=Button(mifr,text="cargar valores",command=guardara)
b1.grid(row=3,column=1,sticky="w")
#estapa de mostrar señal muestreada
b2=Button(mifr,text="Submuestrear las Señales",command=plot)
b2.grid(row=5,column=1,sticky="w")
#etapa de mostrar señal digitalizada
b3=Button(mifr,text="Señales digitalizadas ",command=plot3)
b3.grid(row=6,column=1,sticky="w")
#etapa de ventaneo
b3=Button(mifr,text="DFT ",command=plot4)
b3.grid(row=7,column=1,sticky="w")
#etapa de ventaneo
b3=Button(mifr,text="señales  ",command=valo)
b3.grid(row=8,column=1,sticky="w")

b4=Button(mifr,text="Fasores ",command=fin)
b4.grid(row=9,column=1,sticky="w")




raiz.mainloop()
