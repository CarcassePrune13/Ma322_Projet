# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 08:44:50 2021

@author: carca
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt

def fusee(Y,t):

    D = 4
    a=8*10**3
    g=9.81
    k=0.1
    u=2*10**3
    #print('ddddddd ddddd ',Y[1])
    Yprime = np.zeros(3)
    
    if (Y[1] < 80):
        Y[1] = 80
        D=0
    
    Yprime[0]=D*u/Y[1] -g -k*np.exp(-Y[2]/a)*Y[0]**2/Y[1]
    Yprime[1]=-D
    Yprime[2] = Y[0]
    
    return Yprime

'''
La fusée (odeint)
'''

tf=np.linspace(0,160,100)
N=tf.size


plt.figure(3)
Y0=[0,400,0]
y = odeint(fusee, Y0, tf)

plt.plot(tf,y[:,0],label="Vitesse de la fusée")
plt.xlabel("temps en secondes")
plt.ylabel("Vitesse en m/s")
plt.legend()
plt.grid()
plt.show()

plt.figure(4)
plt.plot(tf,y[:,2],label="Trajectoire")
plt.ylabel("Hauteur en m")
plt.xlabel("temps en secondes")
plt.legend()
plt.grid()
plt.show()