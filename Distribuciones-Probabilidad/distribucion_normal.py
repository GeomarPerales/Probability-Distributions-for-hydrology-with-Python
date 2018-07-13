#DISTRIBUCIONES DE PROBABILIDAD
#DISTRIBUCION NORMAL
import numpy as np
import pandas as pd
from csv import reader
import statistics as stats

data=pd.read_csv('data.csv')
a=data['registro']
Q=data['caudal']

xm=stats.mean(Q) #media
dst=stats.stdev(Q) #desv. standard de la muestra
cv=xm/dst #coef. de variacion

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull.py' A '<title>.py'

Px=[Px22[49],Px22[79],Px22[89],Px22[94],Px22[97],Px22[98],Px33[198],Px44[498],Px55[998]]
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria standard

Q_n=np.zeros((len(VAE)))

#caudales segun Dist. NORMAL para un tiempo de retorno(TR)
TR=[2,5,10,20,50,100,200,500,1000] #tiempo de retorno
for i in range (len(VAE)):
    Q_n[i]=xm+VAE[i]*dst 
