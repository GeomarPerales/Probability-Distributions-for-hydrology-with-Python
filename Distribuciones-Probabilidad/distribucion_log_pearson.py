#hidrologia estadistica
#DISTRIBUCIÓN LOG PEARSON

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
N=len(Q)

lnQ=np.log(Q)
cvQ=dst/xm

xm_lp=stats.mean(lnQ)
ds_lp=stats.stdev(lnQ)
cv_lp=ds_lp/xm_lp

V_LP=[xm,ds_lp,cv_lp]

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull.py' A 'distribucion_<title>'

Px=[Px22[49],Px22[79],Px22[89],Px22[94],Px22[97],Px22[98],Px33[198],Px44[498],Px55[998]]
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria 

#coeficiente asimetria
sg_lp=0
for i in range(len(Q)):
    sg_lp=sg_lp+((np.log(Q[i])-xm_lp)**3)/len(Q)  
    g_lp=len(Q)**2*sg_lp/(len(Q)-1)/(len(Q)-2)/ds_lp**3

cs_lp=sg_lp

gc_lp=cs_lp/(np.sqrt(N*(N-1))/(N-2)*(1+8.5/N))
be_lp=(2/gc_lp)**2 #beta 
sc=ds_lp*np.sqrt(N/(N-1))
al_lp=sc/np.sqrt(be_lp) #alfa
y_lp=xm_lp-al_lp*be_lp

V_LP2=[cs_lp,gc_lp,be_lp,sc,al_lp,y_lp]

#variable t
Q_lpt=np.zeros((len(VAE)))

for i in range(len(VAE)):
    Q_lpt[i]=np.exp(al_lp*be_lp*(1-1/9/be_lp+VAE[i]*(1/9/be_lp)**0.5)**3+y_lp)

#variable k
P111=np.zeros((len(VAE)))
P222=np.zeros((len(VAE)))
kp2=np.zeros((len(VAE)))

for j in range(len(VAE)):
    gc1=gc_lp/6
    P111[j]=VAE[j]+(VAE[j]**2-1)*gc1+(VAE[j]**3-6*VAE[j])*gc1**2/3;
    P222[j]=-(VAE[j]**2-1)*gc1**3+VAE[j]*gc1**4+gc1**5/3;
    kp2[j]=P111[j]+P222[j]

Q_lpk=np.zeros((len(VAE)))

for i in range(len(VAE)):
    Q_lpk[i]=np.exp(xm_lp+kp2[i]*ds_lp)
