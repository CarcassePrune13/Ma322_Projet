# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 08:39:36 2021

@author: carca
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *

C = 10e-6  
#C = 100e-6 #Avec cette valeur, la méthode de Range-Kutta fonctionne
R = 3 
L = 0.5 
e = 10

def RLC_prim_s(Y,t):
    RLC_prim=np.array([Y[1],(1/L)*(e-Y[0]-R*C*Y[1])])
    return (RLC_prim)

def RLC_prim_i(Y,t):
    RLC_prim=np.array([Y[1],-(1/L)*((1/C)*Y[0]+R*Y[1])])
    return (RLC_prim)

def RLC_prim(Y,t):
    RLC_prim=np.array([(1/L)*(e-Y[1]-R*Y[0]),(1/C)*Y[0]])
    return (RLC_prim)

t = np.arange(0,2,0.01)   
Y0 = np.array([0,0])  
N=len(t)
h=0.005
def solv_edo():
    Yode = odeint(RLC_prim, Y0, t)
    return Yode

def RK(f,t,Y0,N,h):
    Yrk=np.zeros((N,Y0.size))
    Yrk[0,:]=Y0.reshape(2)
    for k in range (N-1):
        k1=f(Yrk[k,:],t[k])
        k2=f(Yrk[k,:]+(h/2)*k1,t[k]+(h/2))
        k3=f(Yrk[k,:]+(h/2)*k2,t[k]+(h/2))
        k4=f(Yrk[k,:]+h*k3,t[k]+(h/2))
        Yrk[k+1,:]=Yrk[k,:] + (h/6)*(k1+2*k2+2*k3+k4)
    return (Yrk)

plt.figure(1)
Yrk= RK(RLC_prim,t,Y0,N,h)
Yode=solv_edo()
plt.plot(t,Yrk[:,0],label='Runge Kutta')
plt.plot(t,Yode[:,0],label='Odeint')
plt.title("Système RLC") 
plt.xlabel('temps (s)')
plt.ylabel('s(t)')
plt.grid(True)
plt.legend()
plt.show()