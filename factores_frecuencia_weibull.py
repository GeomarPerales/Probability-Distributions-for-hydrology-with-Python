#hidrologia estadistica
import numpy as np
import math
import astropy
from astropy.table import Table

#TABLA F(Z)
#coeficientes

ao=2.490895
a1=1.466003
a2=-0.024343
a3=0.178257

#tabla de ordenadas
#horizontal
th=np.arange(0,0.1,0.01)

#vertical
tv=np.arange(0,4,0.1)

#matriz ceros
Fz=np.zeros((len(tv),len(th)))

#function
#Fz=(ao+a1*t**2+a2*t**4+a3*t**6)**(-1)
for i in range (len(tv)):
    for j in range(len(th)):
        Fz[i][j]=(ao+a1*(tv[i]+th[j])**2+a2*(tv[i]+th[j])**4+a3*(tv[i]+th[j])**6)**(-1)


#----------------------------------------------------------------------------------------
#ABRAMOWITSY
#calculo de T
#coeficientes
as_a3=0.33267

#matriz zeros
tas=np.zeros((len(tv),len(th)))

#funcion
for i in range (len(tv)):
    for j in range(len(th)):
        tas[i][j]=1/(1+as_a3*(tv[i]+th[j]))
        
#CALCULO DE FZ
#Coeficientes
as_ao=0.43618
as_a1=0.12017
as_a2=0.9373

#nueva vertical
tv1=np.arange(0,3.2,0.1)

Fzas=np.zeros((len(tv1),len(th)))

for i in range (len(tv1)):
    for j in range (len(th)):
        Fzas[i][j]=1-Fz[i][j]*(as_ao*tas[i][j]-as_a1*tas[i][j]**2+as_a2*tas[i][j]**3)

#----------------------------------------------------------------------------------------
#MASTING
#calculo de T
#coeficientes

ma_bo=0.232164

#matriz zeros
tma=np.zeros((len(tv),len(th)))

#funcion
for i in range(len(tv)):
    for j in range(len(th)):
        tma[i][j]=1/(1+ma_bo*(tv[i]+th[j]))
    
#CALCULO DE FZ
#Coeficientes
mb1=0.31938
mb2=-0.35656
mb3=1.78148
mb4=-1.82126
mb5=1.33027

P=np.zeros((len(tv1),len(th)))
Fzma=np.zeros((len(tv1),len(th)))

for i in range (len(tv1)):
    for j in range (len(th)):
        P[i][j]=mb1*tma[i][j]+mb2*tma[i][j]**2+mb3*tma[i][j]**3+mb4*tma[i][j]**4+mb5*tma[i][j]**5
        Fzma[i][j]=1-Fz[i][j]*P[i][j]

#---------------------------------------------------------------------------------------------------
#FACTOR DE FRECUENCIA usando WEIBULL
#---------------------------------------------------------------------------------------------------
#muestra 1:59
N1=np.arange(1,60)

#coeficientes
co=2.515517
c1=0.802853
c2=0.010328

d1=1.432788
d2=0.189269
d3=0.001308

Px1=np.zeros(len(N1))
Px11=np.zeros(len(N1))
TM1=np.zeros(len(N1))
W1=np.zeros(len(N1))
k1=np.zeros(len(N1))

for i in range(len(N1)):
    Px1[i]=(i+1)/(len(N1)+1)
    Px11[i]=(i+1)/(len(N1)+1)
    TM1[i]=1/(1-Px1[i])
    if Px1[i]>=.5:
        Px1[i]=1-Px1[i]
    W1[i]=((math.log(1/Px1[i]**2)))**0.5
    k1[i]=W1[i]-(co+c1*W1[i]+c2*W1[i]**2)/(1+d1*W1[i]+d2*W1[i]**2+d3*W1[i]**3);

FF59=Table([TM1,N1,Px1,W1,k1])



#muestra 1-99
N2=np.arange(1,100)

Px2=np.zeros(len(N2))
Px22=np.zeros(len(N2))
TM2=np.zeros(len(N2))
W2=np.zeros(len(N2))
k2=np.zeros(len(N2))

for i in range(len(N2)):
    Px2[i]=(i+1)/(len(N2)+1)
    Px22[i]=(i+1)/(len(N2)+1)
    TM2[i]=1/(1-Px2[i])
    if Px2[i]>=.5:
        Px2[i]=1-Px2[i]
    W2[i]=((math.log(1/Px2[i]**2)))**0.5
    k2[i]=W2[i]-(co+c1*W2[i]+c2*W2[i]**2)/(1+d1*W2[i]+d2*W2[i]**2+d3*W2[i]**3);

FF99=Table([TM2,N2,Px2,W2,k2])



#muestra 1-199
N3=np.arange(1,200)

Px3=np.zeros(len(N3))
Px33=np.zeros(len(N3))
TM3=np.zeros(len(N3))
W3=np.zeros(len(N3))
k3=np.zeros(len(N3))

for i in range(len(N3)):
    Px3[i]=(i+1)/(len(N3)+1)
    Px33[i]=(i+1)/(len(N3)+1)
    TM3[i]=1/(1-Px3[i])
    if Px3[i]>=.5:
        Px3[i]=1-Px3[i]
    W3[i]=((math.log(1/Px3[i]**2)))**0.5
    k3[i]=W3[i]-(co+c1*W3[i]+c2*W3[i]**2)/(1+d1*W3[i]+d2*W3[i]**2+d3*W3[i]**3);

FF199=Table([TM3,N3,Px3,W3,k3])


            
#muestra 1-499
N4=np.arange(1,500)

Px4=np.zeros(len(N4))
Px44=np.zeros(len(N4))
TM4=np.zeros(len(N4))
W4=np.zeros(len(N4))
k4=np.zeros(len(N4))

for i in range(len(N4)):
    Px4[i]=(i+1)/(len(N4)+1)
    Px44[i]=(i+1)/(len(N4)+1)
    TM4[i]=1/(1-Px4[i])
    if Px4[i]>=.5:
        Px4[i]=1-Px4[i]
    W4[i]=((math.log(1/Px4[i]**2)))**0.5
    k4[i]=W4[i]-(co+c1*W4[i]+c2*W4[i]**2)/(1+d1*W4[i]+d2*W4[i]**2+d3*W4[i]**3);

FF499=Table([TM4,N4,Px4,W4,k4])
            

            
#muestra 1-999
N5=np.arange(1,1000)

Px5=np.zeros(len(N5))
Px55=np.zeros(len(N5))
TM5=np.zeros(len(N5))
W5=np.zeros(len(N5))
k5=np.zeros(len(N5))

for i in range(len(N5)):
    Px5[i]=(i+1)/(len(N5)+1)
    Px55[i]=(i+1)/(len(N5)+1)
    TM5[i]=1/(1-Px5[i])
    if Px5[i]>=.5:
        Px5[i]=1-Px5[i]
    W5[i]=((math.log(1/Px5[i]**2)))**0.5
    k5[i]=W5[i]-(co+c1*W5[i]+c2*W5[i]**2)/(1+d1*W5[i]+d2*W5[i]**2+d3*W5[i]**3);

FF999=Table([TM5,N5,Px5,W5,k5])
            



        
