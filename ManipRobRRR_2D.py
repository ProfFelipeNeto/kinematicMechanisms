"""
Instituto Federal de Educação Ciência e Tecnologia do Rio Grande do Sul
Programa: Cinemática do Manipulador Robotico RRR no plano
Author: Prof. Dr. Felipe Rodrigues de Freitas Neto
Data: 14/09/2023
"""
import numpy as np
from numpy import sin,pi,cos, sqrt
import matplotlib.pyplot as plt

#Parametros
l1, l2, l3= 0.5, 0.3, 0.1 #[m]
thetap_1,thetap_2,thetap_3=5.,10.,25. #rad/s
thetapp_1, thetapp_2, thetapp_3 = 0.0, 0.0, 0.0 #rad/s^2
tempo=0 ; Tempo=[]

dTheta1=0.10*pi/180
#Theta 1
theta1Grad_i=0 #angulo inicial do theta 1 [graus]
theta1Grad_f= -45 #angulo final do theta 1 [graus]
theta1_i= pi*theta1Grad_i/180  #angulo final do theta 1 [rad]
theta1_f= pi*theta1Grad_f/180  #angulo final do theta 1 [rad]
if (theta1Grad_f<0):
    dTheta1=-dTheta1
    thetap_1=-thetap_1
Theta1=np.arange(theta1_i,theta1_f,dTheta1)
TempoTotal_1=np.abs((theta1_f-theta1_i)/thetap_1)
dtempo1=np.abs(dTheta1/thetap_1)

#Theta 2
dTheta2=0.101*pi/180
theta2Grad_i=0 #angulo inicial do theta 2 [graus]
theta2Grad_f= 90 #angulo final do theta 2 [graus]
theta2_i= pi*theta2Grad_i/180
theta2_f= pi*theta2Grad_f/180  #angulo final do theta 2 [rad]
if (theta2Grad_f<0):
    dTheta2=-dTheta2
    thetap_2=-thetap_2
Theta2=np.arange(theta2_i,theta2_f,dTheta2)
TempoTotal_2=np.abs((theta2_f-theta2_i)/thetap_2)
dtempo2=np.abs(dTheta2/thetap_2)

#Theta 3
dTheta3=0.101*pi/180
theta3Grad_i= 90 #angulo inicial do theta 3 [graus]
theta3Grad_f= -90 #angulo final do theta 3 [graus]
theta3_i= pi*theta3Grad_i/180  #angulo inicial do theta 3 [rad]
theta3_f= pi*theta3Grad_f/180  #angulo final do theta 3 [rad]
if (theta3Grad_f<theta3Grad_i):
    dTheta3=-dTheta3 
    thetap_3=-thetap_3
Theta3=np.arange(theta3_i,theta3_f,dTheta3)
TempoTotal_3=np.abs((theta3_f-theta3_i)/thetap_3)
dtempo3=np.abs(dTheta3/thetap_3)

posB=np.empty((0,3),dtype=float) ; posC=np.empty((0,3),dtype=float)
posD=np.empty((0,3),dtype=float)
velB=np.empty((0,4),dtype=float) ; velC=np.empty((0,4),dtype=float)
velD=np.empty((0,4),dtype=float)
acelB=np.empty((0,4),dtype=float) ; acelC=np.empty((0,4),dtype=float)
acelD=np.empty((0,4),dtype=float)

## Funções das Posições
def posicaoB(theta1,tempo):
    rBx= - l1*sin(theta1)
    rBy= l1*cos(theta1)  
    R_ob=np.array([tempo,rBx,rBy])      
    return R_ob

def posicaoC(theta1,theta2,rBx,rBy):
    rCx= rBx + l2*sin(theta1+theta2)
    rCy= rBy -l2*cos(theta1+theta2)  
    R_oc=np.array([theta1_grad,rCx,rCy]) #depois theta1 trocar por tempo    
    return R_oc
     
def posicaoD(theta1,theta2,theta3,rCx,rCy):
    rDx= rCx+l3*cos(theta1+theta2+theta3)
    rDy= rCy+l3*sin(theta1+theta2+theta3)    
    R_od=np.array([theta1_grad,rDx,rDy]) #depois theta1 trocar por tempo    
    return R_od

## Funções das Velocidades    
def velocidadeB(mov,tempo,theta1,thetap_1):
  if (mov==1):
      thetap_1=thetap_1
  else:
      thetap_1=0    
  vBx= -thetap_1*l1*cos(theta1)
  vBy= -thetap_1*l1*sin(theta1)
  
  modVel_B=np.sqrt(vBx**2+vBy**2)
  vB=np.array([tempo,vBx,vBy,modVel_B])
  return vB

def velocidadeC(mov,tempo,theta1,theta2,thetap_1,thetap_2,vBx,vBy) :  
  if(mov==1):
     thetap_1=thetap_1; thetap_2=0
  elif(mov==2):
      thetap_1=0; thetap_2=thetap_2
  elif(mov==3) :
      thetap_1=0; thetap_2=0
          
  vCx=vBx + (thetap_1+thetap_2)*l2*cos(theta1+theta2)
  vCy=vBy + (thetap_1+thetap_2)*l2*sin(theta1+theta2)         
  modVel_C=sqrt(vCx**2+vCy**2)
  vC=np.array([tempo,vCx,vCy,modVel_C])
  return vC

def velocidadeD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,vCx,vCy) :  
  if(mov==1):
     thetap_1=thetap_1; thetap_2=0;thetap_3=0
  elif(mov==2):
      thetap_1=0; thetap_2=thetap_2;thetap_3=0
  elif(mov==3) :
      thetap_1=0; thetap_2=0;thetap_3=thetap_3
        
  vDx= vCx - (thetap_1+thetap_2+thetap_3)*l3*sin(theta1+theta2+theta3)
  vDy= vCy + (thetap_1+thetap_2+thetap_3)*l3*cos(theta1+theta2+theta3)
  modVel_D=sqrt(vDx**2+vDy**2)
  vD=np.array([tempo,vDx,vDy,modVel_D])
  return vD

def aceleracaoB(mov,tempo,theta1,thetap_1,thetapp_1):
   if (mov==1):
      thetap_1=thetap_1
      thetapp_1=thetapp_1
   else:
      thetap_1=0
      thetapp_1=0

   aBx=-thetapp_1*l1*cos(theta1)+(thetap_1)**2*l1*sin(theta1)
   aBy=-thetapp_1*l1*sin(theta1)-(thetap_1)**2*l1*cos(theta1)
   modAcel_B=sqrt(aBx**2+aBy**2)
   aB=np.array([tempo,aBx,aBy,modAcel_B])
   return aB  

def aceleracaoC(mov,tempo,theta1,theta2,thetap_1,thetap_2,thetapp_1, thetapp_2,
                aBx,aBy) :    
    if(mov==1):
       thetap_1=thetap_1; thetap_2=0
       thetapp_1=thetapp_1 ; thetapp_2=0
    elif(mov==2):
        thetap_1=0; thetap_2=thetap_2
        thetapp_1=0 ; thetapp_2=thetapp_2
    elif(mov==3) :
        thetap_1=0; thetap_2=0
        thetapp_1=0 ; thetapp_2=0
         
    aCx=aBx+(thetapp_1+thetapp_2)*l2*cos(theta1+theta2)-(thetap_1+thetap_2)**2*l2*sin(theta1+theta2)
    aCy=aBy+(thetapp_1+thetapp_2)*l2*sin(theta1+theta2)+(thetap_1+thetap_2)**2*l2*cos(theta1+theta2)
    modAcel_C=sqrt(aCx**2+aCy**2)
    aC=np.array([tempo,aCx,aCy,modAcel_C])
    return aC

def aceleracaoD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,thetapp_1,
                thetapp_2,thetapp_3, aCx,aCy) :  
    if(mov==1):
       thetap_1=thetap_1; thetap_2=0; thetap_3=0
       thetapp_1=thetapp_1 ; thetapp_2=0 ;thetapp_3=0
    elif(mov==2):
        thetap_1=0; thetap_2=thetap_2 ; thetap_3=0
        thetapp_1=0 ; thetapp_2=thetapp_2 ;thetapp_3=0
    elif(mov==3) :
        thetap_1=0; thetap_2=0 ; thetap_3=thetap_3
        thetapp_1=0 ; thetapp_2=0 ; thetapp_3=thetapp_3
         
    Stheta=theta1+theta2+theta3
    Sthetap=thetap_1+thetap_2+thetap_3
    Sthetapp=thetapp_1+thetapp_2+thetapp_3
    aDx=aCx -Sthetapp*l3*sin(Stheta) -(Sthetap)**2*l3*cos(Stheta)
    aDy=aCy +Sthetapp*l3*cos(Stheta) -(Sthetap)**2*l3*sin(Stheta)
    modAcel_D=sqrt(aDx**2+aDy**2)
    aD=np.array([tempo,aDx,aDy,modAcel_D])
    return aD

#### Theta 1 ####
for theta1 in Theta1:
  mov=1  
  theta1_grad=(theta1*180/pi)  
  tempo=tempo+dtempo1
  Tempo.append(tempo)
  theta2=theta2_i;   theta3=theta3_i
  # ---------------- POSIÇÃO ----------------
  # Point B
  R_ob=posicaoB(theta1,tempo)    #rBx=R_ob[1], rBy=R_ob[2]
  posB=np.append(posB,[R_ob],axis=0)  
  # Point C
  R_oc=posicaoC(theta1,theta2,R_ob[1],R_ob[2]) #rCx=R_oc[1], rCy=R_oc[2]
  posC=np.append(posC,[R_oc],axis=0)  
  # Point D
  R_od=posicaoD(theta1,theta2,theta3,R_oc[1],R_oc[2])#rDx=R_od[1], rDy=R_od[2]
  posD=np.append(posD,[R_od],axis=0)  
  
  # ---------------- VELOCIDADE ----------------
  # Point B
  vB=velocidadeB(mov,tempo,theta1,thetap_1)
  velB=np.append(velB,[vB],axis=0)    
  # Point C
  vC=velocidadeC(mov,tempo,theta1,theta2,thetap_1,thetap_2,vB[1],vB[2])
  velC=np.append(velC,[vC],axis=0)
  # Point D
  vD=velocidadeD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,vC[1],vC[2])
  velD=np.append(velD,[vD],axis=0)
  
  # ---------------- ACELERAÇÃO ----------------
  # Point B
  aB=aceleracaoB(mov, tempo, theta1, thetap_1, thetapp_1)
  acelB=np.append(acelB,[aB],axis=0)
  # Point C
  aC=aceleracaoC(mov,tempo,theta1,theta2,thetap_1,thetap_2,thetapp_1, thetapp_2,
                  aB[1],aB[2])
  acelC=np.append(acelC,[aC],axis=0)
  # Point D
  aD=aceleracaoD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,thetapp_1,
                  thetapp_2,thetapp_3, aC[1],aC[2])
  acelD=np.append(acelD,[aD],axis=0)
    
#### Theta 2 ####
for theta2 in Theta2:
  mov=2  
  theta2_grad=(theta2*180/pi)  
  theta3=theta3_i
  tempo=tempo+dtempo2
  Tempo.append(tempo)

# ---------------- POSIÇÃO ----------------
  # Point B
  R_ob=posicaoB(theta1,tempo) #rBx=R_ob[1], rBy=R_ob[2]
  posB=np.append(posB,[R_ob],axis=0)  
  # Point C
  R_oc=posicaoC(theta1,theta2,R_ob[1],R_ob[2])  #rCx=R_oc[1], rCy=R_oc[2]
  posC=np.append(posC,[R_oc],axis=0)   
  # Point D
  R_od=posicaoD(theta1,theta2,theta3,R_oc[1],R_oc[2]) #rDx=R_od[1], rDy=R_od[2]
  posD=np.append(posD,[R_od],axis=0)
# ---------------- VELOCITY ----------------
  # Point B
  vB=velocidadeB(mov,tempo,theta1,thetap_1)
  velB=np.append(velB,[vB],axis=0)   
  # Point C
  vC=velocidadeC(mov,tempo,theta1,theta2,thetap_1,thetap_2,vB[1],vB[2])
  velC=np.append(velC,[vC],axis=0)
  # Point D
  vD=velocidadeD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,vC[1],vC[2])
  velD=np.append(velD,[vD],axis=0)
# ---------------- ACELERAÇÃO ----------------
  # Point B
  aB=aceleracaoB(mov, tempo, theta1, thetap_1, thetapp_1)
  acelB=np.append(acelB,[aB],axis=0)
  # Point C
  aC=aceleracaoC(mov,tempo,theta1,theta2,thetap_1,thetap_2,thetapp_1, thetapp_2,
                  aB[1],aB[2])
  acelC=np.append(acelC,[aC],axis=0)
  # Point D
  aD=aceleracaoD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,thetapp_1,
                  thetapp_2,thetapp_3, aC[1],aC[2])
  acelD=np.append(acelD,[aD],axis=0)
  
#### Theta 3 ####  
for theta3 in Theta3:
  mov=3  
  theta3_grad=(theta2*180/pi)  
  tempo=tempo+dtempo3
  Tempo.append(tempo)
  # ---------------- Posição ----------------
  # Point B
  R_ob=posicaoB(theta1,tempo)  
  posB=np.append(posB,[R_ob],axis=0)
  # Point C
  R_oc=posicaoC(theta1,theta2,R_ob[1],R_ob[2])  #rCx=R_oc[1], rCy=R_oc[2]
  posC=np.append(posC,[R_oc],axis=0)    
  # Point D
  R_od=posicaoD(theta1,theta2,theta3,R_oc[1],R_oc[2]) #rDx=R_od[1], rDy=R_od[2]
  posD=np.append(posD,[R_od],axis=0)
 # ---------------- VELOCITY ----------------
  # Point B
  vB=velocidadeB(mov,tempo,theta1,thetap_1)
  velB=np.append(velB,[vB],axis=0)   
  # Point C 
  vC=velocidadeC(mov,tempo,theta1,theta2,thetap_1,thetap_2,vB[1],vB[2])
  velC=np.append(velC,[vC],axis=0)
  # Point D
  vD=velocidadeD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,vC[1],vC[2])
  velD=np.append(velD,[vD],axis=0)  
# ---------------- ACELERAÇÃO ----------------
  # Point B
  aB=aceleracaoB(mov, tempo, theta1, thetap_1, thetapp_1)
  acelB=np.append(acelB,[aB],axis=0)
  # Point C
  aC=aceleracaoC(mov,tempo,theta1,theta2,thetap_1,thetap_2,thetapp_1, thetapp_2,
                  aB[1],aB[2])
  acelC=np.append(acelC,[aC],axis=0)
  # Point D
  aD=aceleracaoD(mov,tempo,theta1,theta2,theta3,thetap_1,thetap_2,thetap_3,thetapp_1,
                  thetapp_2,thetapp_3, aC[1],aC[2])
  acelD=np.append(acelD,[aD],axis=0)  
     
# ####Gráficos####   
#Posição 
plt.figure(figsize=(6,5),dpi=300)
plt.plot(posB[:,1],posB[:,2], linewidth=1.5, color="black",label='Ponto B')
plt.plot(posC[:,1],posC[:,2], linewidth=1.5, color="blue",label='Ponto C')
plt.plot(posD[:,1],posD[:,2], linewidth=1.5, color="red",label='Ponto D')
plt.xlabel('X [m]',fontsize='18') ;plt.ylabel('Y [m]',fontsize='18')
plt.legend(loc='best', fontsize=12) ; plt.xlim(0,0.7); plt.grid()
plt.xticks(fontsize=14); plt.yticks(fontsize=14) ; plt.show()
#Velocidade
plt.figure(figsize=(8,6),dpi=300)
plt.subplot(311)
plt.plot(Tempo,velB[:,3], linewidth=1.5, color="black",label='Mod. vB')
plt.xlim(Tempo[0],Tempo[-1]); plt.legend(loc='best', fontsize=12)
plt.xticks(fontsize=13); plt.yticks(fontsize=13); plt.grid(True)
plt.subplot(312)
plt.plot(Tempo,velC[:,3], linewidth=1.5, color="blue",label='Mod. vC')
plt.ylabel("Velocidade [m/s]",fontsize='18')
plt.xlim(Tempo[0],Tempo[-1]); plt.legend(loc='best', fontsize=12)
plt.xticks(fontsize=13); plt.yticks(fontsize=13); plt.grid(True)
plt.subplot(313)
plt.plot(Tempo,velD[:,3], linewidth=1.5, color="red",label='Mod. vD')
plt.xlabel("Tempo [s]",fontsize='18')
plt.xlim(Tempo[0],Tempo[-1]);plt.ylim(0,1.2*np.max(velD[:,3])); plt.legend(loc='best', fontsize=12); plt.grid(True)
plt.xticks(fontsize=13); plt.yticks(fontsize=13); plt.show()
#Aceleração
plt.figure(figsize=(8,6),dpi=300)
plt.subplot(311)
plt.plot(Tempo,acelB[:,3], linewidth=1.5, color="black",label='Mod. aB')
plt.xlim(Tempo[0],Tempo[-1]); plt.legend(loc='best', fontsize=12); plt.grid(True)
plt.xticks(fontsize=13); plt.yticks(fontsize=13) ; 
plt.subplot(312)
plt.plot(Tempo,acelC[:,3], linewidth=1.5, color="blue",label='Mod. aC')
plt.ylabel("Aceleração [m/s²]",fontsize='18')
plt.xlim(Tempo[0],Tempo[-1]); plt.legend(loc='best', fontsize=12); plt.grid(True)
plt.xticks(fontsize=13); plt.yticks(fontsize=13) ; 
plt.subplot(313)
plt.plot(Tempo,acelD[:,3], linewidth=1.5, color="red",label='Mod. aD')
plt.xlim(Tempo[0],Tempo[-1]); plt.legend(loc='best', fontsize=12); plt.grid(True)
plt.xticks(fontsize=13); plt.yticks(fontsize=13) ; 
plt.xlabel("Tempo [s]",fontsize='18') ;plt.show()