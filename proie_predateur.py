# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 16:09:05 2021

@author: carca
"""

from math import *
from pylab import *
from numpy import *
from scipy.integrate import *

def proies ():
    alph1=3
    L1=list()
    for i in range (0,10):
        if i ==0:
            L1.append(5)
        else:
            L1.append(L1[i-1]*alph1)
    return L1


def predateurs ():
    alph2=-2
    L2=list()
    for i in range (0,10):
        if i ==0:
            L2.append(3)
        else:
            if L2[i-1]<0:
                L2.append(L2[i-1]*-alph2)
            else: 
                L2.append(L2[i-1]*alph2)
    return L2

def Euler(Yprime,t,h,Y0):
    Ye=zeros((N,2))
    Ye[0,:]=Y0
    for i in range(N-1):  
        Ye[i+1,:]=Ye[i,:]+h*Yprime(Ye[i,:],t[i])
    return Ye

def RungeKunta4(Yprime,t,h,Y0):
    YRK=zeros((N,2))
    YRK[0,:]=Y0
    for i in range(N-1):
        Y=YRK[i,:]
        k1=Yprime(Y,t[i])
        k2=Yprime(Y+(h/2)*k1,t[i]+(i/2))
        k3=Yprime(Y+h*(k2/2),t[i]+i/2)
        k4=Yprime(Y+h*k3,t[i]+h)
        YRK[i+1,:]=YRK[i,:]+(h/6)*(k1+2*k2+2*k3+k4)
    return YRK


def proie_predateur(Y,t):
    Yprime = zeros(2)
    alpha1 = 3
    beta1 = 1
    alpha2 = 2
    beta2 = 1

    Yprime[0] = alpha1*Y[0] - beta1*Y[0]*Y[1]
    Yprime[1] = -alpha2*Y[1]+beta2*Y[0]*Y[1]
    return Yprime


L1=proies()
L2=predateurs()


Tp=linspace(0,10,101)
h=10/100
Y0=[5,3]  
N=len(Tp)
Ypp = Euler(proie_predateur,Tp,h,Y0)
YRK=RungeKunta4(proie_predateur,Tp,h,Y0)
Yodep=odeint(proie_predateur,Y0,Tp)
Tpp=linspace(0,10,10)
figure(figsize=(10,5))
plot(Tpp,L1,label="proies",color="blue")
xlabel('années')
ylabel('nombre de spécimen')
legend(loc=1)
grid(True)
title('nombre de proies en fonction du temps sans prédateurs')
show()

figure(figsize=(10,5))
plot(Tpp,L2,label="predateurs",color="red")
xlabel('années')
ylabel('nombre de spécimen')
legend(loc=1)
grid(True)
title('nombre de prédateurs en fonction du temps sans proies')
show()

figure(figsize=(10,5))
plot(Tp,Ypp[:,0], label="proies",color="blue")
plot(Tp,Ypp[:,1], label="prédateurs",color="red")
xlabel("années")
ylabel("nombre d'especes")
grid(True)
legend(loc=1)
title("Euler")
show()

figure(figsize=(10,5))
plot(Tp,YRK[:,0], label="proies",color="blue")
plot(Tp,YRK[:,1], label="prédateurs",color="red")
xlabel("années")
ylabel("nombre d'especes")
grid(True)
legend(loc=1)
title("RungeKutta4")
show()

figure(figsize=(10,5))
plot(Tp,Yodep[:,0], label="proies",color="blue")
plot(Tp,Yodep[:,1], label="prédateurs",color="red")
xlabel("années")
ylabel("nombre d'especes")
grid(True)
legend(loc=1)
title("odeint")
show()

figure(figsize=(10,5))
plot(Yodep[:,0],Yodep[:,1], label="Odeint",color="blue")
plot(YRK[:,0],YRK[:,1], label="Range-Kutta",color="red")
#plot(Ypp[:,0],Ypp[:,1], label="Euler",color="blue")
xlabel("proies")
ylabel("prédateurs")
grid(True)
legend(loc=1)
title("Portrait de phase")
show()
