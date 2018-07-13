#Hidrologia estadistica
#TABLA DE FACTOR DE FRECUENCIA PARA GUMBEL
#PARA PROB. DE 50%,80%,90%, 95% Y 100 CON INC. DE 5 U
import numpy as np
import statistics as stats
import astropy
from astropy.table import Table

#generacion de las series
TM=np.arange(10,105,5)
A=np.zeros((len(TM),100))
n=np.zeros((len(TM)))

for i in range(len(TM)):
    n[i]=len(np.arange(1,(TM[i]+1)))
    for j in range((len(np.arange(TM[i])))):
        A[i][j]=-np.log(-np.log((n[i]+1-(j+1))/(n[i]+1)));

#media y desviación standard de cada serie
M=np.zeros((len(TM)))
DS=np.zeros((len(TM)))
f_TR=np.zeros((len(TM)))

for i in range(len(TM)):
        M[i]=stats.mean(A[i][0:TM[i]])
        DS[i]=stats.pstdev(A[i][0:TM[i]])

tabFF=Table([TM,M,DS]) #tabla media y desv. standard

TR=[2,5,10,20,50,100,200,500,1000] #tiempo de retorno

for i in range(len(TR)):
    f_TR[i]=-np.log(-np.log((TR[i]-1)/TR[i]))

Yt=np.zeros((len(TM),len(f_TR)))

for i in range(len(TM)):
    for j in range(len(f_TR)):
        Yt[i][j]=-(M[i]-f_TR[j])/DS[i]
