# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 08:38:23 2021

@author: carca
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def NoyauT(x,t): 
    y=2 
    return y

def FT(t):
    y=cos(pi*t/2) - 8/pi 
    return y

def Mat(f,K,a,b,N):
    h=(b-a)/(N-1)
    x = np.linspace(a,b,N)
    t = np.linspace(a,b,N)
    F = np.zeros(N)
    A = np.zeros((N,N))
    for i in range(N):
        F[i]=f(t[i])
        for j in range(N):
            A[i,j] = 2*K(x[i],t[j])
    
    A[:,0] = A[:,0]/2
    A[:,-1] = A[:,-1]/2
    F = F.T
    M = np.diag(np.ones(N)) -(h/2)*A
    
    return (t,F, M)

#Exploitation du code equation intégrale

t,F,M = Mat(FT,NoyauT,-1,1,10)
U = np.linalg.solve(M,F)
exact=[]
for k in range (len(t)):
    ex = cos(pi*t[k]/2)
    exact.append(ex)
plt.figure(0)
plt.plot(t,U,label="Approchée")
plt.plot(t,exact,label="Exact")
plt.grid()
plt.legend(loc = 'best') 
plt.show()
err=np.linalg.norm((U-exact),2)

print("Erreur : ", err)

#Equation de Love

def f(t):
    y=1
    return y

def Noyau(x,t):
    y=(1/pi)*(1/(1+((x-t)**2)))
    return y

t,F,M = Mat(f,Noyau,-1,1,10)
U = np.linalg.solve(M,F)
print("Valeur numérique de U : ",U)
plt.figure(1)
plt.plot(t,U,label="Approchée")

plt.grid()
plt.legend(loc = 'best') 
plt.show()