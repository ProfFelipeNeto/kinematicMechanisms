"""
Instituto Federal do Rio Grande do Sul
Programa: Cinemática de um Guindaste 
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 04/10/2023
"""
import numpy as np
from numpy import cos,sin,sqrt,pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

tfinal,dt= 3,0.005 #[s]
T=np.arange(0,tfinal+dt,dt)
#Parametros
l=40.0 # [m]

omega1, omega2= 0.5, 1.5#[rad/s]
alpha1, alpha2 = 0, 0 #[rad/s^2]

theta1=0*(pi/180) ; theta1_f=45*(pi/180); theta1_old=theta1
theta2=0*(pi/180) ; theta2_f=60*(pi/180); theta2_old=theta2

posB=np.empty((0,7),dtype=float) ; velB=np.empty((0,9),dtype=float)
acelB=np.empty((0,7),dtype=float)
#Posição
def posicaoB(t,theta1,theta2):
    rOBx= l*cos(theta1)*sin(theta2)
    rOBy= l*sin(theta1)*sin(theta2)
    rOBz= l*cos(theta2)   
    modB=sqrt(rOBx**2+rOBy**2+rOBz**2)
    rB=np.array([t,rOBx,rOBy,rOBz,modB,(theta1*180/pi),(theta2*180/pi) ])    
    return rB
#Velocidade
def velocidadeB(t,theta1,theta2,omega1,omega2):
    vBx= l*(omega2*cos(theta1)*cos(theta2) - omega1*sin(theta1)*sin(theta2))
    vBy= l*(omega2*sin(theta1)*cos(theta2) + omega1*cos(theta1)*sin(theta2))
    vBz= -l*omega2*sin(theta2)   
    modvB=sqrt(vBx**2+vBy**2+vBz**2)
    vB=np.array([t,vBx,vBy,vBz,modvB,(theta1*180/pi),(theta2*180/pi),omega1,omega2])     
    return vB
#Aceleração
def aceleracaoB(t,theta1,theta2,omega1,omega2,alpha1,alpha2):
    aTanx= l*(alpha2*cos(theta1)*cos(theta2) -alpha1*sin(theta1)*sin(theta2))
    aTany= l*(alpha2*sin(theta1)*cos(theta2) +alpha1*cos(theta1)*sin(theta2))
    aTanz= -l*alpha2*sin(theta2)    
    aNormx= -l*((omega1**2+omega2**2)*cos(theta1)*sin(theta2) + omega1*omega2*sin(theta1)*cos(theta2))  
    aNormy= -l*((omega1**2+omega2**2)*sin(theta1)*sin(theta2) - omega1*omega2*cos(theta1)*cos(theta2)) 
    aNormz= -l*(omega2)**2*cos(theta2)
    
    aBx= aTanx + aNormx ;     aBy= aTany + aNormy ;   aBz= aTanz + aNormz
    modaB=sqrt(aBx**2+aBy**2+aBz**2)
    aB=np.array([t,aBx,aBy,aBz,modaB,(theta1*180/pi),(theta2*180/pi)])      
    return aB

for t in T:       
        
    #Posição do ponto B
    rB=posicaoB(t,theta1,theta2)   
    posB=np.append(posB,[rB],axis=0)
    #Velocidade do ponto B
    vB=velocidadeB(t,theta1,theta2,omega1,omega2)
    velB=np.append(velB,[vB],axis=0)
    #Aceleração do ponto B
    aB=aceleracaoB(t,theta1,theta2,omega1,omega2,alpha1,alpha2)
    acelB=np.append(acelB,[aB],axis=0)
        
    #Theta 1
    dTheta_1=omega1*dt 
    if (omega1<0 and (abs(theta1)>(theta1_f-abs(dTheta_1)) and abs(theta1)<(theta1_f+abs(dTheta_1)))):
        theta1=theta1_f        
    elif(omega1<0 and (abs(theta1)<(theta1_f-abs(dTheta_1)) or abs(theta1)>(theta1_f+abs(dTheta_1)))):
        theta1=theta1+dTheta_1        
    elif(omega1>0 and (theta1>(theta1_f-abs(dTheta_1)) and theta1<(theta1_f+abs(dTheta_1)))): #ok
        theta1=theta1_f  
    elif(omega1>0 and (theta1<(theta1_f-abs(dTheta_1)) or theta1>(theta1_f+abs(dTheta_1)))): #ok
        theta1=theta1+dTheta_1    
    #Theta 2
    dTheta_2=omega2*dt
    if (omega2<0 and (abs(theta2)>(theta2_f-abs(dTheta_2)) and abs(theta2)<(theta2_f+abs(dTheta_2)))):
      theta2=theta2_f        
    elif(omega2<0 and (abs(theta2)<(theta2_f-abs(dTheta_2)) or abs(theta2)>(theta2_f+abs(dTheta_2)))):
      theta2=theta2+dTheta_2        
    elif(omega2>0 and (theta2>(theta2_f-abs(dTheta_2)) and theta2<(theta2_f+abs(dTheta_2)))): #ok
      theta2=theta2_f  
    elif(omega2>0 and (theta2<(theta2_f-abs(dTheta_2)) or theta2>(theta2_f+abs(dTheta_2)))): #ok
      theta2=theta2+dTheta_2     
    #omega 1         
    if(theta1==theta1_old):  omega1=0 
    else:   omega1=omega1+alpha1*dt        
    #omega 2         
    if(theta2==theta2_old):  omega2=0
    else:  omega2=omega2+alpha2*dt  
    theta1_old=theta1; theta2_old=theta2
    #alpha 1         
    if(theta1==theta1_old):  alpha1=0 
    else:   alpha1=alpha1        
    #alpha 2         
    if(theta2==theta2_old):  alpha2=0
    else:  alpha2=alpha2    
    theta1_old=theta1; theta2_old=theta2

#Gráfico da Posição    ;ax.set_zlim(10,30) 
fig=plt.figure(figsize=(8,6),dpi=300)    
ax=fig.add_subplot(111,projection='3d')
ax.plot(posB[:,1].flatten(),posB[:,2].flatten(),posB[:,3].flatten(), label='Ponto B',color='black',linewidth=1.5)
ax.set_xlabel('X [m]',fontsize=14) ; ax.set_ylabel('Y [m]',fontsize=14) ; ax.set_zlabel('Z [m]',fontsize=14)
ax.view_init(elev=30, azim=70) ; plt.show(block=True)
#Gráfico da Velocidade
plt.figure(figsize=(6,4),dpi=300)
plt.plot(velB[:,0],velB[:,4],label='mod. vB',color='black',linewidth=1.5) ; plt.grid(True)
plt.xlim(velB[0,0],velB[-1,0]); plt.xlabel("Tempo [s]",fontsize='18'); plt.ylabel("Velocidade [m/s]",fontsize='18')
plt.xticks(fontsize=13); plt.yticks(fontsize=13); plt.legend(loc='best', fontsize=12) ;plt.show()
#Gráfico da Aceleração
plt.figure(figsize=(6,4),dpi=300)
plt.plot(acelB[:,0],acelB[:,4],label='mod. aB',color='black',linewidth=1.5) ; plt.grid(True)
plt.xlim(acelB[0,0],acelB[-1,0]); plt.xlabel("Tempo [s]",fontsize='18'); plt.ylabel("Aceleração [m/s²]",fontsize='18')
plt.xticks(fontsize=13); plt.yticks(fontsize=13); plt.legend(loc='best', fontsize=12) ;plt.show()    