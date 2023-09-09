"""
Instituto Federal de Educação Ciência e Tecnologia do Rio Grande do Sul
Programa: Cinemática do Mecanismo Biela Manivela
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 09/09/2023
"""

import numpy as np
from numpy import sin,pi,cos
import matplotlib.pyplot as plt

#Parametros

L1 = 0.1 #[m]
L2 = 0.4 #[m]

theta1=0.0
dtheta1=0.1*pi/180
n=1200 #rpm
omega1= n*2*pi/60 #rad/s
alpha1=0.0 # rad/s²

#dt=dTheta/omega1
#T=np.arange(0,10,dt) #[s]

Y=[]
Theta1=np.arange(theta1,4*pi,dtheta1)

Theta1_grad=[]
Theta2=[]

Vcy=[]
Acy=[]

Omega1=[]
Omega2=[]
Alpha1=[]
Alpha2=[]


for theta1 in Theta1:
    
    
    #Posição
    theta2 = np.arcsin(L1*cos(theta1)/L2)
    y=L1*sin(theta1)+L2*cos(theta2)
   
    
    Y.append(y)
    Theta2.append(theta2)
    
    theta1_grad=theta1*180/pi
    Theta1_grad.append(theta1_grad)
    
    #Velocidade
    omega2=-omega1*L1*sin(theta1)/(L2*cos(theta2))
    vcy=omega1*L1*cos(theta1) - omega2*L2*sin(theta2)
    
    Vcy.append(vcy)
    Omega1.append(omega1)
    Omega2.append(omega2)
    
    #Aceleração
    alpha2=(-alpha1*L1*sin(theta1) - omega1**2*L1*cos(theta1) + omega2**2*L2*sin(theta2)) / (L2*cos(theta2))
    acy= alpha1*L1*cos(theta1)-omega1**2*L1*sin(theta1)-omega2**2*L2*sin(theta2)-alpha2*L2*sin(theta2)
    
    Acy.append(acy)
    Alpha1.append(alpha1)
    Alpha2.append(alpha2)
    
    

####Gráficos####

#Posição
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta1_grad,Y, color='red',linewidth=1.5)
plt.xlim(0,720)
plt.ylabel('Y [m]', fontsize=14)
plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.grid(True)

#Velocidade
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta1_grad,Vcy, color='red',linewidth=1.5)
plt.xlim(0,720)
plt.ylabel('Vcy [m/s]', fontsize=14)
plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.grid(True)

#Aceleração
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta1_grad,Acy, color='red',linewidth=1.5)
plt.ylabel('Acy [m/s²]', fontsize=14)
plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.xlim(0,720)
plt.grid(True)

plt.show()