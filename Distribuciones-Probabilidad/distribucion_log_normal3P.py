#hidrologia estadistica
#LOG NORMAL 3P
import numpy as np
import pandas as pd
from csv import reader
import statistics as stats
import math

data=pd.read_csv('data.csv')
a=data['registro']
Q=data['caudal']

xm=stats.mean(Q)
dst=stats.stdev(Q)
cv=dst/xm

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull.py' A 'distribucion_<title>'

VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]]

#coeficiente asimetria
sg=0

for i in range(len(Q)):
    sg=sg+((Q[i]-xm)**3)/len(Q)

g=len(Q)**2*sg/(len(Q)-1)/(len(Q)-2)/dst**3
cs=g

W=(-g+math.sqrt(g**2+4))*0.5
Z2=(1-W**(2/3))/W**(1/3)
dy2=(np.log(Z2**2+1))**(1/2)
uy2=np.log(dst/Z2)-0.5*np.log(Z2**2+1)
xo=xm-dst/Z2

V_LN3P=[cs,W,Z2,dy2,uy2,xo]

Q_ln3t=np.zeros((len(VAE)))
k_cv3=np.zeros((len(VAE)))
Q_ln3k=np.zeros((len(VAE)))

#usando t - hallamos Q
for i in range(len(VAE)):
    Q_ln3t[i]=xo+math.exp(uy2+VAE[i]*dy2) #caudales segun Dist. LOG PEARSON 3P

#usando k - hallamos Q
for i in range(len(VAE)):
    k_cv3[i]=(np.exp((np.log(1+Z2**2))**0.5*VAE[i]-0.5*(np.log(1+Z2**2)))-1)/Z2

for i in range(len(VAE)):
    Q_ln3k[i]=xm+k_cv3[i]*dst #caudales segun Dist. LOG PEARSON 3P

