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
        
        
