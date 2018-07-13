#hidrologia estadistica
#DISTRIBUCION PEARSON

import numpy as np
import statistics as stats
import pandas as pd
from csv import reader

data=pd.read_csv('data.csv')
a=data['registro']
Q=data['caudal']

xm=stats.mean(Q)
dst=stats.stdev(Q)
cv=dst/xm

sg=0
#coeficiente asimetria
for i in range(len(Q)):
    sg=sg+((Q[i]-xm)**3)/len(Q)

g=len(Q)**2*sg/(len(Q)-1)/(len(Q)-2)/dst**3
cs=g

#variables a usar
#media xm
#desv. std d_std
#coef. variacion cv
#coef. simetria g
N=len(Q)
gc=cs/(np.sqrt(N*(N-1))/(N-2)*(1+8.5/N))
be=(2/gc)**2 #beta 
alf=dst/(be)**0.5  #alfa
y=xm-dst*(be)**0.5;

VP=[gc,be,alf,y]

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull.py' A 'distribucion_<title>.py'
Px=[Px22[49],Px22[79],Px22[89],Px22[94],Px22[97],Px22[98],Px33[198],Px44[498],Px55[998]]
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria
#mediante t
Q_pt=np.zeros((len(Px)))

for i in range(len(Px)):
    Q_pt[i]=alf*be*(1-1/9/be+VAE[i]*math.sqrt(1/9/be))**3+y

#mediante k
P11=np.zeros((len(VAE)))
P22=np.zeros((len(VAE)))
kp1=np.zeros((len(VAE)))

for j in range(len(VAE)):
    gc1=gc/6
    P11[j]=VAE[j]+(VAE[j]**2-1)*gc1+(VAE[j]**3-6*VAE[j])*gc1**2/3
    P22[j]=-(VAE[j]**2-1)*gc1**3+VAE[j]*gc1**4+gc1**5/3
    kp1[j]=P11[j]+P22[j]

Q_pk=np.zeros((len(VAE)))

for i in range(len(VAE)):
    Q_pk[i]=xm+kp1[i]*dst
