#Hidrologia estadistica
#DISTRIBUCION GUMBEL

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

al=1.2825/dst #alfa
u=xm-0.45*dst #u

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull' A 'distribucion_<title>.py'
Px=[Px22[49],Px22[79],Px22[89],Px22[94],Px22[97],Px22[98],Px33[198],Px44[498],Px55[998]]

Ytg=np.zeros((len(Px)))
Q_gt=np.zeros((len(Px)))

#mediante t
for i in range(len(Px)):
    Ytg[i]=-np.log(-np.log(Px[i]))
    Q_gt[i]=Ytg[i]/al+u

#mediante k
A59=np.zeros((len(Q)))
    
for i in range(len(Q)):
        A59[i]=-np.log(-np.log((len(Q)-i)/(len(Q)+1)))
 
xmg=stats.mean(A59)
d_stg=stats.pstdev(A59)

VG=[xmg,d_stg]

k_g=np.zeros((len(Ytg)))
Q_gk=np.zeros((len(Ytg)))

for i in range(len(Ytg)):
        k_g[i]=(Ytg[i]-xmg)/d_stg
        Q_gk[i]=xm+k_g[i]*dst
