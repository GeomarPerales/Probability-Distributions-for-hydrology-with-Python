#Hidrologia estadistica
#TABLA DE FRECUENCIA PEARSON
#prob .... y coe;f de sesgo de 0.0 a 2.0 con inc de 0.1

import numpy as np

csp=np.arange(0,2.1,0.1)

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull' A 'distribucion_<title>.py'
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria standard

gcp=np.zeros((len(csp)))
P1=np.zeros((len(csp),len(VAE)))
P2=np.zeros((len(csp),len(VAE)))
k_p=np.zeros((len(csp),len(VAE)))

for i in range(len(csp)):
    for j in range(len(VAE)):
        gcp[i]=csp[i]/6;
        P1[i][j]=VAE[j]+(VAE[j]**2-1)*gcp[i]+(VAE[j]**3-6*VAE[j])*gcp[i]**2/3
        P2[i][j]=-(VAE[j]**2-1)*gcp[i]**3+VAE[j]*gcp[i]**4+gcp[i]**5/3
        k_p[i][j]=P1[i][j]-P2[i][j]
