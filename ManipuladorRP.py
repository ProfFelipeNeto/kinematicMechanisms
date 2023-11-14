"""
Instituto Federal do Rio Grande do Sul
Programa: Cinemática do Manipulador RP
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 14/11/2023
"""
import numpy as np 
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
from matplotlib.pyplot import grid,legend,xlabel,ylabel,xlim, xticks, yticks

t=0 ; tFinal=5 ; dt=0.1 #[s]
T=np.arange(t,tFinal,dt)
omega= 2*pi/16 #[rad/s]
alpha=0 #[rad/s^2]
thetaInicial=0 ; theta=thetaInicial + omega*t #[rad]
L= 0.4 + 0.1*t #[m]
Lp=0.1 #[m/s]
Lpp=0 # [m/s^2]

posB=np.empty((0,3),dtype=float) ; VelB=np.empty((0,4),dtype=float)
AcelB=np.empty((0,4),dtype=float)

for t in T:
    
  #Posição
    #Condição para o avanço do braço
    if(L>=0.7):   L=0.7
    else:    L= 0.4 + 0.1*t
    #Condição para o movimento rotativo
    if(theta>pi):  theta=pi
    elif(theta<0): theta=0
    else: theta=thetaInicial + omega*t
   
    rBx=L*cos(theta) ;     rBy=L*sin(theta)
    posb=np.array([t,rBx,rBy]); posB=np.append(posB,[posb],axis=0)
  #Velocidade
    xVb=-omega*L*sin(theta) + Lp*cos(theta) 
    yVb= omega*L*cos(theta) + Lp*sin(theta)  
    modVb=np.sqrt(xVb**2+yVb**2)
    velB=np.array([t,xVb,yVb,modVb]); VelB=np.append(VelB, [velB],axis=0)
  #Aceleração
    xAb= -alpha*L*sin(theta)-omega**2*L*cos(theta)-2*omega*Lp*sin(theta)+Lpp*cos(theta)
    yAb=  alpha*L*cos(theta)-omega**2*L*sin(theta)+2*omega*Lp*cos(theta)+Lpp*sin(theta)
    modAb=np.sqrt(xAb**2+yAb**2)
    acelB=np.array([t,xAb,yAb,modAb])
    AcelB=np.append(AcelB,[acelB],axis=0)

#Gráficos
#Posição
plt.figure(figsize=(6,6),dpi=300)
plt.plot(posB[:,1],posB[:,2],label='Ponto B',color="red",linewidth=1.5)
grid(True); xlabel('X [m]',fontsize=18);ylabel('Y [m]',fontsize=18)
legend(loc='best',fontsize=14);xticks(fontsize=12);yticks(fontsize=12)
plt.figure(figsize=(8,3),dpi=300)
plt.plot(posB[:,0],posB[:,1],label='xB',color="black",linewidth=1.5)
plt.plot(posB[:,0],posB[:,2],label='yB',color="blue",linewidth=1.5)
grid(True);xlabel('Tempo [s]',fontsize=18);ylabel('Posição [m]',fontsize=18)
xlim(posB[0,0],posB[-1,0]); xticks(fontsize=14); yticks(fontsize=14)
legend(loc='best',fontsize=14); plt.show()
#Velocidade
plt.figure(figsize=(8,6),dpi=300)
plt.subplot(211)
plt.plot(VelB[:,0],VelB[:,1],label='xVb',color="black",linewidth=1.5)
plt.plot(VelB[:,0],VelB[:,2],label='yVb',color="blue",linewidth=1.5)
grid(True);ylabel('Velocidade [m/s]',fontsize=18); yticks(fontsize=14)
xlim(posB[0,0],posB[-1,0]); xticks(fontsize=14);legend(loc='best',fontsize=14)
plt.subplot(212)
plt.plot(VelB[:,0],VelB[:,3],label='mod Vb',color="red",linewidth=1.5)
grid(True);xlabel('Tempo [s]',fontsize=18);ylabel('Velocidade [m/s]',fontsize=18)
xlim(VelB[0,0],VelB[-1,0]); xticks(fontsize=14); yticks(fontsize=14)
legend(loc='best',fontsize=14); plt.show()
#Aceleração
plt.figure(figsize=(8,6),dpi=300)
plt.subplot(211)
plt.plot(AcelB[:,0],AcelB[:,1],label='xAb',color="black",linewidth=1.5)
plt.plot(AcelB[:,0],AcelB[:,2],label='yAb',color="blue",linewidth=1.5)
grid(True);ylabel('Aceleração [m/s²]',fontsize=16); yticks(fontsize=14)
xlim(AcelB[0,0],AcelB[-1,0]); xticks(fontsize=14);legend(loc='best',fontsize=14)
plt.subplot(212)
plt.plot(AcelB[:,0],AcelB[:,3],label='mod Ab',color="red",linewidth=1.5)
grid(True);xlabel('Tempo [s]',fontsize=18);ylabel('Aceleração [m/s²]',fontsize=16)
xlim(AcelB[0,0],AcelB[-1,0]); xticks(fontsize=14); yticks(fontsize=14)
legend(loc='best',fontsize=14); plt.show()