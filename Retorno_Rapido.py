"""
Instituto Federal do Rio Grande do Sul
Programa: Cinemática do Mecanismo de Retorno Rápido
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 13/08/2023
"""

import numpy as np
from numpy import sin,pi,cos, tan
import matplotlib.pyplot as plt

#Parametros
L1, L4, L5 = 0.1, 0.4, 0.2 #[m]

n=300 #rpm
omega1= n*2*pi/60 #rad/s
alpha1=0.0 # rad/s²
dTheta=0.101*pi/180
Theta=np.arange(0,4*pi,dTheta) ; Theta_grad=[]

# Posição
Phi=[] ; L2=[]
Roc=np.empty((0,2),dtype=float) ; Rod=np.empty((0,2),dtype=float)
#Velocidade
Phip=[] #Derivada primeira de phi [rad/s]
L2p = [] #Derivada primeira de L2 [m/s]
Vc=np.empty((0,3),dtype=float) ; Vd=np.empty((0,3),dtype=float)
#Aceleração
Phipp=[] #Derivada segunda de phi [rad/s²]
L2pp = [] #Derivada segunda de L2 [m/s²]
Ac=np.empty((0,3),dtype=float) ; Ad=np.empty((0,3),dtype=float)

for theta in Theta:    
    theta_grad=(theta*180/pi)
    
 #**************POSIÇÃO **************
 #Ponto C   
    Rocx=L1*cos(theta) 
    Rocy=L1*sin(theta)    
    r_oc=np.array([Rocx,Rocy]);  Roc=np.append(Roc, [r_oc],axis=0)    
    phi = np.arctan((L1*sin(theta)+L4)/(L1*cos(theta)))
    l2= L1*cos(theta)/cos(phi)    
 #Ponto D   
    l3=(L4+L5)/sin(phi)    
    Rodx=l3*cos(phi)     
    Rody=L5
    r_od=np.array([Rodx,Rody]); Rod=np.append(Rod,[r_od],axis=0)    
    Phi.append(phi) ;     L2.append(l2)

#**************VELOCIDADE **************    
 #Velocidade C   
    phip=(omega1*L1*cos(theta-phi))/(l2) ; Phip.append(phi)
    l2p=omega1*L1*cos(theta)/sin(phi) - phip*l2/tan(phi) ; L2p.append(l2p)       
    Vcx=-omega1*L1*sin(theta)
    Vcy= omega1*L1*cos(theta)
    mod_vc=np.sqrt(Vcx**2+Vcy**2)    
    vc=np.array([Vcx,Vcy,mod_vc]) ;  Vc=np.append(Vc,[vc],axis=0)    
 #Velocidade D
    l3p= -phip*l3/tan(phi)
    Vdx= -phip*l3*sin(phi) + l3p*cos(phi)
    Vdy=0.0
    mod_Vd=np.sqrt(Vdx**2+Vdy**2)    
    vd=np.array([Vdx,Vdy,mod_Vd]) ;  Vd=np.append(Vd,[vd],axis=0)

#**************ACELERAÇÃO **************   
#Aceleração C   
    b1= -omega1**2*L1*cos(theta) + phip**2*l2*cos(phi) + 2*phip*l2p*sin(phi)
    b2= -omega1**2*L1*sin(theta) + phip**2*l2*sin(phi) + 2*phip*l2p*sin(phi)
    phipp= -(b1*sin(phi) - b2*cos(phi))/l2
    l2pp= b2*sin(phi) + b1*cos(phi)
    Acx= -omega1**2*L1*cos(theta) 
    Acy= -omega1**2*L1*sin(theta)
    mod_Ac=np.sqrt(Acx**2+Acy**2)    
    ac=np.array([Acx,Acy,mod_Ac]) ;    Ac=np.append(Ac,[ac],axis=0)
#Aceleração D
    l3pp= -phipp*l3/tan(phi) +phip**2*l3 -2*phip*l3p/tan(phi)
    Adx= -phipp*l3*sin(phi) -phip**2*l3*cos(phi) -2*phip*l3p*sin(phi) + l3pp*cos(phi)
    Ady= 0.0    
    mod_Ad=np.sqrt(Adx**2+Ady**2)    
    ad=np.array([Adx,Ady,mod_Ad]) ;    Ad=np.append(Ad,[ad],axis=0)
    Theta_grad.append(theta_grad)
   
############ GRÁFICOS DE  POSIÇÃO ############
plt.figure(figsize=(8,5),dpi=300)
#Ponto C
plt.subplot(211)
plt.plot(Theta_grad,Roc[:,0], linewidth=1.5, color="black",label='Xc')
plt.plot(Theta_grad,Rod[:,0], linewidth=1.5, color="red",label='Xd')
plt.legend(loc='best') ; plt.xlim(0,720) ; plt.grid(True)
plt.ylabel('Posição X [m]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.subplot(212)
plt.plot(Theta_grad,Roc[:,1], linewidth=1.5, color="black",label='Yc')
plt.plot(Theta_grad,Rod[:,1], linewidth=1.5, color="red",label='Yd')
plt.legend(loc='best') ; plt.xlim(0,720) ; plt.grid(True)
plt.ylabel('Posição Y [m]', fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.xlabel('Theta [$^o$] ', fontsize=14) ; plt.show()

############ GRÁFICOS DE  VELOCIDADE ############
#Ponto C
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plt.plot(Theta_grad,Vc[:,0], linewidth=1.5, color="black",label='Vcx')
plt.plot(Theta_grad,Vc[:,1], linewidth=1.5, color="red",label='Vcy')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Velocidade [m/s]', fontsize=14)
plt.subplot(212)
plt.plot(Theta_grad,Vc[:,2], linewidth=1.5, color="blue",label='Vc (Módulo)')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Velocidade [m/s]', fontsize=14)
plt.xlabel('Theta [$^o$] ', fontsize=14) ; plt.show()
#Ponto D
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plt.plot(Theta_grad,Vd[:,0], linewidth=1.5, color="black",label='Vdx')
plt.plot(Theta_grad,Vd[:,1], linewidth=1.5, color="red",label='Vdy')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Velocidade [m/s]', fontsize=14)
plt.subplot(212)
plt.plot(Theta_grad,Vd[:,2], linewidth=1.5, color="blue",label='Vd (Módulo)')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Velocidade [m/s]', fontsize=14)
plt.xlabel('Theta [$^o$] ', fontsize=14) ; plt.show()

############ GRÁFICOS DE  ACELERAÇÃO ############
#Ponto C
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plt.plot(Theta_grad,Ac[:,0], linewidth=1.5, color="black",label='Acx')
plt.plot(Theta_grad,Ac[:,1], linewidth=1.5, color="red",label='Acy')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Aceleração [m/s²]', fontsize=14)
plt.subplot(212)
plt.plot(Theta_grad,Ac[:,2], linewidth=1.5, color="blue",label='Ac (Módulo)')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Aceleração [m/s²]', fontsize=14)
plt.xlabel('Theta [$^o$] ', fontsize=14) ; plt.show()
#Ponto D
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plt.plot(Theta_grad,Ad[:,0], linewidth=1.5, color="black",label='Adx')
plt.plot(Theta_grad,Ad[:,1], linewidth=1.5, color="red",label='Ady')
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.ylabel('Aceleração [m/s²]', fontsize=14)
plt.subplot(212)
plt.plot(Theta_grad,Ad[:,2], linewidth=1.5, color="blue",label='Ad (Módulo)')
plt.ylabel('Aceleração [m/s²]', fontsize=14)
plt.xlim(0,720) ; plt.legend(loc='best') ; plt.grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True)
plt.xlabel('Theta [$^o$] ', fontsize=14) ; plt.show()
