# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 09:50:53 2021

@author: carca
"""
import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
R=5#ohm
L=50*10**(-3)#Henry
Ke=0.2#V.rad-1.s-1
Kc=0.1#N.m.A-1
Fm=0.01#N.m.rad-1.s
Jm=0.05#kh.m2

l=list()

def Ut(t): #Fonction échelon de tension
    if (t<10 or t>50):
        return 0
    else:
        return 5
    
t=np.arange(0,80,0.01)

for i in t :
    l.append(Ut(i))
    
def moteurCC(Y,t):
    return np.array([(1/Jm)*(Kc*Y[1] - Fm*Y[0]),(1/L)*(Ut(t)-R*Y[1]-Ke*Y[0])])

Yode = odeint(moteurCC,[0,0],t)

plt.plot(t,Yode[:,0],'y',label="Odeint")
plt.title("Evolution de la vitesse angulaire")
plt.xlabel("t en secondes")
plt.ylabel("w(t) en rad/s")
plt.legend()
plt.grid()
plt.show()

plt.plot(t,Kc*Yode[:,1],'g',label="Odeint")
plt.title("Evolution du couple moteur")
plt.xlabel("t en secondes")
plt.ylabel("N.m")
plt.legend()
plt.grid()
plt.show()


plt.plot(t,l,'b',label="Tension")
plt.title("Evolution de la tension u(t)")
plt.xlabel("t en secondes")
plt.ylabel("u(t) en V")
plt.legend()
plt.grid()
plt.show()
