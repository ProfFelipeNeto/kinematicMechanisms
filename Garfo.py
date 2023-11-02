"""
Instituto Federal do Rio Grande do Sul
Programa: Cinemática do Mecanismo Garfo Escocês
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 28/10/2023
"""
import numpy as np
from numpy import sin,pi,cos
import matplotlib.pyplot as plt

#Parametros
L1, L2, h = 0.1,0.4, 0.04 #[m] ; 
theta=0.0
dtheta=0.1*pi/180
n=1200 #rpm ; 
omega= n*2*pi/60 #rad/s
alpha=0.0 # rad/s²
Theta=np.arange(theta,4*pi,dtheta)

xC=[]; Vcx=[] ; Acx=[] ;Theta_grad=[]; yBC=[]

for theta in Theta:
        
    #Posição
    xc = L1*cos(theta)+L2 ;     xC.append(xc)
    y_bc=h-L1*sin(theta)     ;  yBC.append(y_bc)  
    theta_grad=theta*180/pi
    Theta_grad.append(theta_grad)    
    #Velocidade
    yp_bc= -omega*L1*cos(theta) 
    vcx= -omega*L1*sin(theta) 
    Vcx.append(vcx)    
    #Aceleração
    acx= -alpha*L1*sin(theta)-omega**2*L1*cos(theta)    
    ypp_bc= -alpha*L1*cos(theta)+omega**2*L1*sin(theta)        
    Acx.append(acx)

#####Gráficos####
#Posição
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta_grad,xC, color='red',linewidth=1.5)
plt.ylabel('Xc [m]', fontsize=14) ; plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True) ;plt.xlim(0,720) ; 
plt.grid(True)
#Velocidade
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta_grad,Vcx, color='red',linewidth=1.5)
plt.ylabel('Vc [m/s]', fontsize=14) ;plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.xlim(0,720) ; plt.grid(True)
#Aceleração
plt.figure(figsize=(8,3),dpi=300)
plt.plot(Theta_grad,Acx, color='red',linewidth=1.5)
plt.ylabel('Acx [m/s²]', fontsize=14) ; plt.xlabel('Theta [$^0$]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.xlim(0,720) ; plt.grid(True) ; plt.show()