"""
Instituto Federal do Rio Grande do Sul
Programa: Cinemática do Mecanismo de Quatro Barras
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 09/11/2023
"""
import numpy as np
from numpy import pi,cos,sin
import matplotlib.pyplot as plt
from matplotlib.pyplot import grid,xlim,legend,plot,xlabel,ylabel
import random
from scipy.optimize import fsolve

posB=np.empty((0,3),dtype=float); posC=np.empty((0,3),dtype=float)
velB=np.empty((0,4),dtype=float); velC=np.empty((0,4),dtype=float)
acelB=np.empty((0,4),dtype=float);acelC=np.empty((0,4),dtype=float)

theta_1=0.01*pi/180 ; dTheta_1=0.1*pi/180 ; thetaFinal_1=4*pi
Theta1=np.arange(theta_1,thetaFinal_1,dTheta_1)

#Parametros
l1,l2,l3,l4=0.1,0.4,0.3,0.3
n=1200 #rpm
omega1= n*2*pi/60 #rad/s
alpha1=0 #rad/s^2
#Equação de Grashof
Max=l2 ; Min=l1
if ((Max+Min)>(l3+l4)):
  print("Os comprimentos dos elos não satisfaz a condição de Grashof!")  

#Equação não linear
def NonLinEqua(variables,theta_1):
    theta_2,theta_3=variables
    f1 = l2*cos(theta_2)-l3*cos(theta_3)+l1*cos(theta_1)-l4
    f2 = l2*sin(theta_2)-l3*sin(theta_3)+l1*sin(theta_1)    
    return [f1, f2]

#Estimativas iniciais para theta_2 e theta_3
theta_2=random.random() ; theta_3=random.random()

for theta_1 in Theta1:
#Posição do Ponto B      
  xB=l1*cos(theta_1) ;   yB=l1*sin(theta_1)
  thetaGraus_1=theta_1*180/pi
  pB=np.array(([thetaGraus_1,xB,yB]))
  posB=np.append(posB,([pB]),axis=0)
#Posição do Ponto C
  estimInicial = np.array([theta_2, theta_3])  
  theta_2,theta_3 = fsolve(NonLinEqua,estimInicial,args=(theta_1),xtol=1e-08)
  xC1=xB+l2*cos(theta_2) ;   yC1=yB+l2*sin(theta_2)  
  pC=np.array([thetaGraus_1,xC1,yC1])
  posC=np.append(posC, ([pC]),axis=0)
#Velocidade do Ponto B
  xVb= -omega1*l1*sin(theta_1)
  yVb=  omega1*l1*cos(theta_1) 
  modVB= np.sqrt(xVb**2+yVb**2)
  velb=np.array([thetaGraus_1,xVb,yVb,modVB])
  velB=np.append(velB,[velb],axis=0)
#Velocidade do Ponto C
  omega2=(omega1*l1*l3*sin(theta_3-theta_1))/(l2*l3*sin(theta_2-theta_3))
  omega3=(omega1*l1*l2*sin(theta_2-theta_1))/(l2*l3*sin(theta_2-theta_3))  
  xVc=xVb-omega2*l2*sin(theta_2)
  yVc=yVb+omega2*l2*cos(theta_2)
  modVc=np.sqrt(xVc**2+yVc**2)  
  velc=np.array([thetaGraus_1,xVc,yVc,modVc])
  velC=np.append(velC,[velc],axis=0)  
#Aceleração do Ponto B
  xAb= -omega1**2*l1*cos(theta_1)  
  yAb= -omega1**2*l1*sin(theta_1)  
  modAB=np.sqrt(xAb**2+yAb**2)
  acelb=np.array([thetaGraus_1,xAb,yAb,modAB])
  acelB=np.append(acelB,[acelb],axis=0)  
#Aceleração do Ponto C
  detA=l2*l3*sin(theta_2-theta_3)
  b1= omega1**2*l1*cos(theta_1)+omega2**2*l2*cos(theta_2)-omega3**2*l3*cos(theta_3) 
  b2= omega1**2*l1*sin(theta_1)+omega2**2*l2*sin(theta_2)-omega3**2*l3*sin(theta_3) 
  alpha2=-l3*(b1*cos(theta_3)+b2*sin(theta_3))/detA
  alpha3=-l2*(b1*cos(theta_2)+b2*sin(theta_2))/detA
  xAc= xAb -alpha2*l2*sin(theta_2)-omega2**2*l2*cos(theta_2)
  yAc= yAb +alpha2*l2*cos(theta_2)-omega2**2*l2*sin(theta_2)
  modAC=np.sqrt(xAc**2+yAc**2)
  acelc=np.array([thetaGraus_1,xAc,yAc,modAC])
  acelC=np.append(acelC,[acelc],axis=0)  
    
#-------- Gráficos de Posição --------
plt.figure(figsize=(7,6),dpi=300)
plot(posB[:,1],posB[:,2],linewidth=1.5, color="black",label='posB')
plot(posC[:,1],posC[:,2],linewidth=1.5, color="red",label='posC')
ylabel("Y [m]",fontsize=16); xlabel("X [m]",fontsize=16)
legend(loc='best', fontsize=14);grid(True); plt.xticks(fontsize=14)
plt.yticks(fontsize=14); plt.show()

plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plot(posB[:,0],posB[:,1],linewidth=1.5, color="black",label='xB')
plot(posB[:,0],posB[:,2],linewidth=1.5, color="blue",label='yB')
ylabel("Posição [m]",fontsize=14); legend(loc='best')
plt.xticks(np.arange(0,720,step=20),minor=True);xlim(0,720);grid(True)
plt.subplot(212)
plot(posC[:,0],posC[:,1],linewidth=1.5, color="black",label='xC')
plot(posC[:,0],posC[:,2],linewidth=1.5, color="blue",label='yC')
ylabel("Posição [m]",fontsize=14); xlabel("Theta [$^o$]",fontsize=14)
plt.xticks(np.arange(0,720,step=20),minor=True)
xlim(0,720); legend(loc='best');grid(True); plt.show()
#-------- Gráficos de Velocidade --------
#Ponto B
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plot(velB[:,0],velB[:,1],linewidth=1.5, color="black",label='vBx')
plot(velB[:,0],velB[:,2],linewidth=1.5, color="blue",label='vBy')
ylabel("Velocidade [m/s]",fontsize=14); legend(loc='best')
plt.xticks(np.arange(0,720,step=20),minor=True); xlim(0,720); grid(True)
plt.subplot(212)
plot(velB[:,0],velB[:,3],linewidth=1.5, color="red",label='mod. vB')
ylabel("Velocidade [m/s]",fontsize=14); legend(loc='best')
xlabel("Theta [$^o$]",fontsize=14);xlim(0,720); grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True); plt.show()
#Ponto C
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plot(velC[:,0],velC[:,1],linewidth=1.5, color="black",label='vCx')
plot(velC[:,0],velC[:,2],linewidth=1.5, color="blue",label='vCy')
ylabel("Velocidade [m/s]",fontsize=14); legend(loc='best')
plt.xticks(np.arange(0,720,step=20),minor=True); xlim(0,720); grid(True)
plt.subplot(212)
plot(velC[:,0],velC[:,3],linewidth=1.5, color="red",label='mod. vC')
ylabel("Velocidade [m/s]",fontsize=14); legend(loc='best')
xlabel("Theta [$^o$]",fontsize=14); xlim(0,720); grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True);plt.show()
#-------- Gráficos de Aceleração --------
#Ponto B
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plot(acelB[:,0],acelB[:,1],linewidth=1.5, color="black",label='aBx')
plot(acelB[:,0],acelB[:,2],linewidth=1.5, color="blue",label='aBy')
ylabel("Aceleração [m/s²]",fontsize=14); legend(loc='best')
plt.xticks(np.arange(0,720,step=20),minor=True);xlim(0,720); grid(True)
plt.subplot(212)
plot(acelB[:,0],acelB[:,3],linewidth=1.5, color="red",label='mod. aB')
ylabel("Aceleração [m/s²]",fontsize=14); legend(loc='best')
xlabel("Theta [$^o$]",fontsize=14);xlim(0,720); grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True); plt.show()
#Ponto C
plt.figure(figsize=(8,5),dpi=300)
plt.subplot(211)
plot(acelC[:,0],acelC[:,1],linewidth=1.5, color="black",label='aCx')
plot(acelC[:,0],acelC[:,2],linewidth=1.5, color="blue",label='aCy')
ylabel("Aceleração [m/s²]",fontsize=14); legend(loc='best')
plt.xticks(np.arange(0,720,step=20),minor=True);xlim(0,720); grid(True)
plt.subplot(212)
plot(acelC[:,0],acelC[:,3],linewidth=1.5, color="red",label='mod. aC')
ylabel("Aceleração [m/s²]",fontsize=14); legend(loc='best')
xlabel("Theta [$^o$]",fontsize=14);xlim(0,720); grid(True)
plt.xticks(np.arange(0,720,step=20),minor=True); plt.show()