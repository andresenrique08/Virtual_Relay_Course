# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 06:57:05 2020

@author: anavr
"""

import matplotlib.pyplot as plt
import pyComtrade
import numpy as np
from scipy import signal, fftpack
from tkinter import *
from tkinter import filedialog
from inspect import *

# app
#funcion para cargar

global comtradeObj
comtradeObj = pyComtrade.ComtradeRecord()

#------------------ FUNCIONES DEL PROGRAMA -------------------------------
# Para abrir archivo cfg
def abrir_cfg():
    
    global ar
    
    ar=filedialog.askopenfilename(title="Abrir cfg")#,filetypes=(("Archivos COMTRADE .cfg y .dat",".")))  
    
    #fs=
    messagebox.showinfo("Se cargó el archivo .cfg de la siguiente dirección:\n", ar)
    
# Para abrir archivo dat
def abrir_dat():
    
    global dat
    
    dat=filedialog.askopenfilename(title="Abrir cfg")#,filetypes=(("Archivos COMTRADE .cfg y .dat",".")))  
    
    #fs=
    messagebox.showinfo("Se cargó el archivo .dat de la siguiente dirección:\n", dat)
    
# Carga las entradas del programa
def cargar():    
    global e1,e2,e3,e4,e5, ratio_V, ratio_C
    e1=en1.get()
    e2=en2.get()
    e3=en3.get()
    e4=en4.get()
    e5=en5.get()
    
    if e1 == "Seleccionar":
        messagebox.showerror("Error","Debe seleccionar la frecuencia de muestro del relé")
    else:
        e1 = int(e1)
    
    if e2=="" and e3=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación de voltaje")
    elif e2=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación del lado primario de voltaje")
    elif e3=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación del lado secundario de voltaje")
    else:
        e2 = int(e2)
        e3 = int(e3)
        
    if e4=="" and e5=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación de corriente")
    elif e4=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación del lado primario de corriente")
    elif e5=="":
        messagebox.showerror("Error","Debe ingresar la relación de transformación del lado secundario de corriente")
    else:
        e4 = int(e4)
        e5 = int(e5)
        
    if side.get() == 0:
        messagebox.showerror("Error","Debe seleccionar el lado del transformador que desea visualizar")
    elif side.get()==1:
        ratio_V = e2/e3
        ratio_C = e4/e5
    elif side.get()==2:
        ratio_V = 1
        ratio_C = 1
        
    DFP()
    
# Funciones para cada etapa
# ETAPA DE SUBMUESTREO
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

# ETAPA DE DIGITALIZACION
def quantizer(data, quantizing_bits):
    # data is the vector with the signal
    # quantizing_bits is the number of bits for the converter
    # Quantizer - S&H and ADC
    quantizing_levels   = 2 ** quantizing_bits
    quantizing_step     = (np.max(data)-np.min(data)) / quantizing_levels
    quantizing_signal   = np.round (data / quantizing_step) * quantizing_step;
    
    return quantizing_signal

#ETAPA DE DFT
def DFT(time, data, fk, fs_user_cycle):
    # time is the vector of time
    # data is the vector with the signal
    # fk is the frequency of the system
    # fs_user_cycle is the sample rate given by user
    
    N=np.int(fs_user_cycle)
    N_tot = len(data)-N
    Xc = [0]*N_tot
    Xs = [0]*N_tot
    t = [0]*N_tot
    
    for i in np.arange(N_tot):
        xi=data[i:i+N]
        t[i]=time[i]
        Xc_sum = 0
        Xs_sum = 0
        for k in np.arange(N):
            Xc_temp=xi[k]*np.cos(2*np.pi*k/(N))
            Xc_sum=Xc_sum+Xc_temp
            Xs_temp=xi[k]*np.sin(2*np.pi*k/(N))
            Xs_sum=Xs_sum+Xs_temp
            
        Xc[i]= 2/(N*np.sqrt(2))*Xc_sum
        Xs[i]= 2/(N*np.sqrt(2))* Xs_sum
        
    return t, Xc, Xs
        
# Realiza todo el proceso de DSP
def DFP():
    # Definición de variables globales para resultados
    global time, voltages, currents, time_sub, V_sub, C_sub, dig_V_sub, dig_C_sub, t, X_V, X_C, Y_V, Y_C, Xc_V, Xs_V, Xc_C, Xs_C, fs_user_cycle
    
    comtradeObj.read(ar,dat)
    N = comtradeObj['endsamp'][-1]
    #sampling_freq=comtradeObj['samp'][-1]
    fs_comtrade=comtradeObj['samp'][-1]
    fk = comtradeObj['line_freq']
    time = comtradeObj.get_timestamps()
    voltages = np.empty(([len(time),3]))
    currents = np.empty(([len(time),3]))
    # Reading voltaje and currents
    voltages[:,0] = comtradeObj['A'][16]['values']
    voltages[:,1] = comtradeObj['A'][17]['values']
    voltages[:,2] = comtradeObj['A'][18]['values']
    
    currents[:,0] = comtradeObj['A'][0]['values']
    currents[:,1] = comtradeObj['A'][1]['values']
    currents[:,2] = comtradeObj['A'][2]['values']
    for i in np.arange(6):
        if i<3:
            # voltages[:,i] = comtradeObj['A'][i]['values']
            for j in np.arange(len(voltages[:,i])):
                voltages[j,i]= voltages[j,i] * ratio_V
        else:
            # currents[:,i-3] = comtradeObj['A'][i]['values']
            for j in np.arange(len(currents[:,i-3])):
                currents[j,i-3]= currents[j,i-3] * ratio_C
            
    # Submuestreo
    fs_user_cycle = e1 # ESTA DATO VIENE DE LA INTERFAZ
    N_tot = np.int(N*fk/fs_comtrade)*fs_user_cycle
    V_sub = np.empty(([N_tot,3]))
    C_sub = np.empty(([N_tot,3]))
    time_sub = np.empty(([N_tot,6]))
    for i in np.arange(6):
        if i<3:
            time_sub[:,i], V_sub[:,i] = subsampling(time,voltages[:,i],fs_comtrade,fk,fs_user_cycle)
        else:
            time_sub[:,i], C_sub[:,i-3] = subsampling(time,currents[:,i-3],fs_comtrade,fk,fs_user_cycle)
    
    # Digitalización
    quantizing_bits_V = 12 # Valor típico: 12 (Voltaje)
    quantizing_bits_C = 16 # Valor típico: 16 (Corriente)
    dig_V_sub = np.empty(([N_tot,3]))
    dig_C_sub = np.empty(([N_tot,3]))
    for i in np.arange(6):
        if i<3:
            dig_V_sub[:,i] = quantizer(V_sub[:,i], quantizing_bits_V)
        else:
            dig_C_sub[:,i-3] = quantizer(C_sub[:,i-3], quantizing_bits_C)
    
    # DFT
    N_tot_DTF = np.int(N_tot-fs_user_cycle)
    Xc_V = np.empty(([N_tot_DTF,3]))
    Xs_V = np.empty(([N_tot_DTF,3]))
    Xc_C = np.empty(([N_tot_DTF,3]))
    Xs_C = np.empty(([N_tot_DTF,3]))
    X_V = np.empty(([N_tot_DTF,3]))
    Y_V = np.empty(([N_tot_DTF,3]))
    X_C = np.empty(([N_tot_DTF,3]))
    Y_C = np.empty(([N_tot_DTF,3]))
    t = np.empty(([N_tot_DTF,6]))
    for i in np.arange(6):
        if i<3:
            t[:,i], Xc_V[:,i], Xs_V[:,i] = DFT(time_sub[:,i], dig_V_sub[:,i], fk, fs_user_cycle)
            X_V[:,i] = np.sqrt(np.power(Xc_V[:,i],2)+np.power(Xs_V[:,i],2))
            ajus = np.pi
            if Xc_V[-1,i]>0 and Xs_V[-1,i]<0:
                ajus = 2*np.pi
            elif Xc_V[-1,i]>0 and Xs_V[-1,i]>0:
                ajus = 0
            Y_V[:,i] = (np.arctan(Xs_V[:,i]/Xc_V[:,i])+ajus)*180/np.pi
        else:
            t[:,i], Xc_C[:,i-3], Xs_C[:,i-3] = DFT(time_sub[:,i], dig_C_sub[:,i-3], fk, fs_user_cycle)
            X_C[:,i-3] = np.sqrt(np.power(Xc_C[:,i-3],2)+np.power(Xs_C[:,i-3],2))
            ajus = np.pi
            if Xc_C[-1,i-3]>0 and Xs_C[-1,i-3]<0:
                ajus = 2*np.pi
            elif Xc_C[-1,i-3]>0 and Xs_C[-1,i-3]>0:
                ajus = 0
            Y_C[:,i-3] = (np.arctan(Xs_C[:,i-3]/Xc_C[:,i-3])+ajus)*180/np.pi
    
# ------------------------- Funciones para los Botones --------------------
def seniales_COMTRADE():    
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
     
def submuestreo_boton():
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
    
def digitalizacion_boton():
    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('Digitalización de la señal', y=0.92, fontsize=16)
    
    # Plot Voltages
    axarr[0,0].plot(time_sub[:,0][0:fs_user_cycle*3], V_sub[:,0][0:fs_user_cycle*3], 'b-', label='Phase A')
    axarr[0,0].plot( time_sub[:,0][0:fs_user_cycle*3], dig_V_sub[:,0][0:fs_user_cycle*3], 'c-', label='Phase A digital')
    axarr[1,0].plot(time_sub[:,1][0:fs_user_cycle*3], V_sub[:,1][0:fs_user_cycle*3], 'r-', label='Phase B')
    axarr[1,0].plot( time_sub[:,1][0:fs_user_cycle*3], dig_V_sub[:,1][0:fs_user_cycle*3], 'm-', label='Phase B digital')
    axarr[2,0].plot(time_sub[:,2][0:fs_user_cycle*3], V_sub[:,2][0:fs_user_cycle*3], 'g-', label='Phase C')
    axarr[2,0].plot( time_sub[:,2][0:fs_user_cycle*3], dig_V_sub[:,2][0:fs_user_cycle*3], 'y-', label='Phase C digital')
    
    # Plot Currents
    axarr[0,1].plot(time_sub[:,3][0:fs_user_cycle*3], C_sub[:,0][0:fs_user_cycle*3], 'b-', label='Phase A')
    axarr[0,1].plot( time_sub[:,3][0:fs_user_cycle*3], dig_C_sub[:,0][0:fs_user_cycle*3], 'c-', label='Phase A digital')
    axarr[1,1].plot(time_sub[:,4][0:fs_user_cycle*3], C_sub[:,1][0:fs_user_cycle*3], 'r-', label='Phase B')
    axarr[1,1].plot( time_sub[:,4][0:fs_user_cycle*3], dig_C_sub[:,1][0:fs_user_cycle*3], 'm-', label='Phase B digital')
    axarr[2,1].plot(time_sub[:,5][0:fs_user_cycle*3], C_sub[:,2][0:fs_user_cycle*3], 'g-', label='Phase C')
    axarr[2,1].plot( time_sub[:,5][0:fs_user_cycle*3], dig_C_sub[:,2][0:fs_user_cycle*3], 'y-', label='Phase C digital')
    
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
        
def DFT_boton():
    # PLOTING -----------------------------------------------------------------
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('DFT En Magnitud', y=0.92, fontsize=16)
    
    # Plot Voltages
    #axarr[0,0].plot(time_sub[:,0], V_sub[:,0], 'b-', label='Phase A')
    axarr[0,0].plot( t[:,0], X_V[:,0], 'c-', label='Phase A FFT(mag)')
    #axarr[1,0].plot(time_sub[:,1], V_sub[:,1], 'r-', label='Phase B')
    axarr[1,0].plot( t[:,1], X_V[:,1], 'm-', label='Phase B FFT(mag)')
    #axarr[2,0].plot(time_sub[:,2], V_sub[:,2], label='Phase C')
    axarr[2,0].plot( t[:,2], X_V[:,2], 'y-', label='Phase C FFT(mag)')
    
    # Plot Currents
    #axarr[0,1].plot(time_sub[:,3], C_sub[:,0], 'b-', label='Phase A')
    axarr[0,1].plot( t[:,3], X_C[:,0], 'c-', label='Phase A FFT(mag)')
    #axarr[1,1].plot(time_sub[:,4], C_sub[:,1], 'r-', label='Phase B')
    axarr[1,1].plot( t[:,4], X_C[:,1], 'm-', label='Phase B FFT(mag)')
    #axarr[2,1].plot(time_sub[:,5], C_sub[:,2], 'g-', label='Phase C')
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
    
    # Ploting angle
    f, axarr = plt.subplots(3, 2, figsize =(16, 10))
    f.suptitle('DFT En Fase', y=0.92, fontsize=16)
    
    # Plot Voltages
    #axarr[0,0].plot(time_sub[:,0], V_sub[:,0], 'b-', label='Phase A')
    axarr[0,0].plot( t[:,0], Y_V[:,0], 'c-', label='Phase A FFT[ang(rad)]')
    #axarr[1,0].plot(time_sub[:,1], V_sub[:,1], 'r-', label='Phase B')
    axarr[1,0].plot( t[:,1], Y_V[:,1], 'm-', label='Phase B FFT[ang(rad)]')
    #axarr[2,0].plot(time_sub[:,2], V_sub[:,2], label='Phase C')
    axarr[2,0].plot( t[:,2], Y_V[:,2], 'y-', label='Phase C FFT[ang(rad)]')
    
    # Plot Currents
    #axarr[0,1].plot(time_sub[:,3], C_sub[:,0], 'b-', label='Phase A')
    axarr[0,1].plot( t[:,3], Y_C[:,0], 'c-', label='Phase A FFT[ang(rad)]')
    #axarr[1,1].plot(time_sub[:,4], C_sub[:,1], 'r-', label='Phase B')
    axarr[1,1].plot( t[:,4], Y_C[:,1], 'm-', label='Phase B FFT[ang(rad)]')
    #axarr[2,1].plot(time_sub[:,5], C_sub[:,2], 'g-', label='Phase C')
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
    
def fasores_boton():
    # Creando la figura
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
    
    # Mostrando las fases en la interfaz
    label_fas0.grid(row=3,column=0)
    label_fas1.config(text=['Voltaje fase A:', "{:.2f}".format(np.double(X_V[-1:,0])), 'V', "{:.2f}".format(np.double(Y_V[-1:,0])),'°'])
    label_fas1.grid(row=4,column=0)
    label_fas2.config(text=['Voltaje fase B:', "{:.2f}".format(np.double(X_V[-1:,1])), 'V', "{:.2f}".format(np.double(Y_V[-1:,1])),'°'])
    label_fas2.grid(row=5,column=0)
    label_fas3.config(text=['Voltaje fase C:', "{:.2f}".format(np.double(X_V[-1:,2])), 'V', "{:.2f}".format(np.double(Y_V[-1:,1])),'°'])
    label_fas3.grid(row=6,column=0)
    label_fas4.config(text=['Corriente fase A:', "{:.2f}".format(np.double(X_C[-1:,0])), 'A', "{:.2f}".format(np.double(Y_C[-1:,0])),'°'])
    label_fas4.grid(row=7,column=0)
    label_fas5.config(text=['Corriente fase B:', "{:.2f}".format(np.double(X_C[-1:,1])), 'A', "{:.2f}".format(np.double(Y_C[-1:,1])),'°'])
    label_fas5.grid(row=8,column=0)
    label_fas6.config(text=['Corriente fase C:', "{:.2f}".format(np.double(X_C[-1:,2])), 'A', "{:.2f}".format(np.double(Y_C[-1:,2])),'°'])
    label_fas6.grid(row=9,column=0)
#------------------------------------------------------------------
# MAIN
raiz=Tk()
raiz.title("Modulo DSP para relés")
    #raiz.geometry("650x380")
x_frame = 720
y_frame = 500
raiz.geometry("{width}x{height}".format(width=x_frame, height=y_frame)) 
raiz.resizable(True, True)   
# mifr= Frame()
# mifr.pack(side=TOP, fill=BOTH, expand=Y)
# #mifr.config()
# mifr.config(width="650",height="380")
# mifr.config(cursor="star")
bmenu=Menu(raiz)
raiz.config(menu=bmenu)

#-------------------------------- FRAME 1 -------------------------------
frame1 = LabelFrame(raiz, text="Entradas",height= y_frame,width =x_frame,padx=5, labelanchor=N)
frame1.config(cursor="star")
frame1.pack(expand = 'no', fill = 'both') 

x_frame11 = 350
frame11 = LabelFrame(frame1, text="Selección del archivo COMTRADE",fg="red",height= 100,width =x_frame11,padx=15)
frame11.grid(column=0, row=0, padx=10, pady=10)
label11 = Label(frame11, text = 'Instrucciones: Se deben seleccionar los archivos .cfg y .dat')
label11.place(x = 0, y = 5) 

BotonAbrircfg = Button(frame11,text="Abrir .cfg",command=abrir_cfg)
BotonAbrircfg.place(x = x_frame11/6, y = 30) 
BotonAbrirdat = Button(frame11,text="Abrir .dat",command=abrir_dat)
BotonAbrirdat.place(x = 3*x_frame11/6, y = 30) 

x_frame12 = 250
frame12 = LabelFrame(frame1, text="Selección de parámetros de entrada",fg="red",height= 200,width =x_frame12,padx=15)
frame12.grid(column=1, row=0, columnspan=1, rowspan=2, padx=10, pady=10)

t1=Label(frame12,text="Seleccione la frecuencia \nde muestro del relé",fg="green")
t1.grid(row=0,column=0,sticky="w",pady="20")
en1=StringVar()
d=OptionMenu(frame12, en1, "4","8","16","32","64")
d.grid(row=0,column=1,sticky="w",padx="20",pady="20")
en1.set("Seleccionar")


frame121 = LabelFrame(frame12, text="Relación de transformadores de instrumentación",fg="Blue",height= 200,width =x_frame12,padx=15)
frame121.grid(column=0, row=1, columnspan=2, rowspan=1, padx=10, pady=10)

# Relación de transformadores de instrumentación voltaje le falta dividir entre prim y secundario
t2=Label(frame121,text="Relacion de Voltaje")
t2.grid(row=0,column=0,sticky="w")
en2=Entry(frame121, width=7)
en2.grid(row=0,column=1)
t3=Label(frame121,text=":")
t3.grid(row=0,column=2)
en3=Entry(frame121, width=7)
en3.grid(row=0,column=3)

# Relación de transformadores de instrumentación corriente le falta dividir entre prim y secundario
t21=Label(frame121,text="Relacion de Corriente")
t21.grid(row=1,column=0,sticky="w")
en4=Entry(frame121, width=7)
en4.grid(row=1,column=1)
t31=Label(frame121,text=":")
t31.grid(row=1,column=2)
en5=Entry(frame121, width=7)
en5.grid(row=1,column=3)

b1=Button(frame12,text="cargar valores",command=cargar)
b1.grid(row=2,column=0,columnspan=2, rowspan=1)

x_frame13 = 350
frame13 = LabelFrame(frame1, text="Visualización en el transformador",fg="red",height= 100,width =x_frame13,padx=15)
frame13.grid(column=0, row=1, padx=10, pady=0)
label13 = Label(frame13, text = 'Seleccione el lado que desea ver la señal del transformador')
label13.place(x = 0, y = 5) 

side = IntVar()
rad_trafo1 = Radiobutton(frame13,text="Primario", variable = side, value=1)
rad_trafo1.place(x = x_frame13/6, y = 30) 
rad_trafo2 = Radiobutton(frame13,text="Secundario", variable = side, value=2)
rad_trafo2.place(x = 3*x_frame13/6, y = 30) 

#-------------------------------- FRAME 2 -------------------------------
frame2 = LabelFrame(raiz, text="Resultados del procesamiento",height= y_frame,width =x_frame,padx=5, labelanchor=N)
frame2.config(cursor="star")
frame2.pack(expand = 'no', fill = 'both') 

x_frame21 = 350
frame21 = LabelFrame(frame2, text="Detalle por etapa",fg="red",height= 100,width =x_frame21,padx=15)
frame21.grid(column=0, row=0, padx=10, pady=10)
label21 = Label(frame21, text = 'Seleccione la etapa que desea detallar')
label21.grid(row=0,column=0,columnspan=3)

bi=Button(frame21,text="Señales de entrada",command=seniales_COMTRADE)#carga los archivos
bi.grid(row=1,column=0)
b2=Button(frame21,text="Submuestreo",command=submuestreo_boton)
b2.grid(row=1,column=1)
#etapa de mostrar señal digitalizada
b3=Button(frame21,text="Digitalizacion",command=digitalizacion_boton)
b3.grid(row=1,column=2)
#etapa de ventaneo
b3=Button(frame21,text="DFT",command=DFT_boton)
b3.grid(row=2,column=1)



x_frame22 = 250
frame22 = LabelFrame(frame2, text="Señal Procesada",fg="red",height= 100,width =x_frame22,padx=15)
frame22.grid(column=1, row=0, padx=10, pady=10)
label22 = Label(frame22, text = 'En este módulo se muestran las señales resultantes')
label22.grid(row=0,column=0)

b3=Button(frame22,text="Fasores",command=fasores_boton)
b3.grid(row=1,column=0)

label_fas0 = Label(frame22, text = 'Los fasores resultantes en rms son:')
label_fas0.grid_forget()
label_fas1 = Label(frame22, text = 'Aqui se mostrara el fasor 1',fg="blue")
label_fas1.grid_forget()
label_fas2 = Label(frame22, text = 'Aqui se mostrara el fasor 2',fg="red")
label_fas2.grid_forget()
label_fas3 = Label(frame22, text = 'Aqui se mostrara el fasor 3',fg="green")
label_fas3.grid_forget()
label_fas4 = Label(frame22, text = 'Aqui se mostrara el fasor 1',fg="blue")
label_fas4.grid_forget()
label_fas5 = Label(frame22, text = 'Aqui se mostrara el fasor 2',fg="red")
label_fas5.grid_forget()
label_fas6 = Label(frame22, text = 'Aqui se mostrara el fasor 3',fg="green")
label_fas6.grid_forget()



barchiv=Menu(bmenu,tearoff=0)
bhelp=Menu(bmenu,tearoff=0)
barchiv.add_command(label="Abrir archivo ")#carga los archivos
barchiv.add_command(label="Guardar")
bmenu.add_cascade(label="Archivo",menu=barchiv)
help

bhelp.add_command(label="Descripción")
bmenu.add_cascade(label="Ayuda",menu=bhelp)



raiz.mainloop()


