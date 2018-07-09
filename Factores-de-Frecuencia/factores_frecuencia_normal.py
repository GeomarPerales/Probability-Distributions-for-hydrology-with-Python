#DISTRIBUCIONES DE PROBABILIDAD
#TABLA DE FACTOR DE FRECUENCIA para prob. 50%, 80%, 90% y 100% y
#coeficiente de variación (CV) que varian de 0.05 hasta 1.00.

import numpy as np
import math

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull.py' A '<title>.py'
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria standard


cvn=np.arange(0.05,1.05,0.05)
kn=np.zeros((len(cvn),len(VAE)))

for i in range (len(cvn)):
    for j in range(len(VAE)):
        kn[i][j]=(math.exp((math.log(1+cvn[i]**2))**0.5*VAE[j]-0.5*(math.log(1+cvn[i]**2)))-1)/cvn[i]
