#hidrologia estadistica

import numpy as np
importh math

#LOG NORMAL 2P
lnQ=np.log(Q)
cvQ=dst/xm

uy1=0.5*math.log(xm**2/(1+cvQ**2))
dy1=(math.log(1+cvQ**2))**0.5

V_DN2P=[cvQ,uy1,dy1]

#valores obtenidos de la tabla de frecuencias de weibull
#NOTA: ACOPLAR 'factores_frecuencia_weibull' A 'distribucion_<title>'
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]]

Q_ln2t=np.zeros((len(VAE)))
cv2p=np.zeros((len(VAE)))
Q_ln2k=np.zeros((len(VAE)))

#sando VAE t
for i in range(len(VAE)):
    Q_ln2t[i]=math.exp(uy1+VAE[i]*dy1); #caudales segun Dist. LOG PEARSON 2P

#usando k
for i in range(len(VAE)):
    cv2p[i]=(np.exp((np.log(1+cvQ**2))**0.5*VAE[i]-0.5*(np.log(1+cvQ**2)))-1)/cvQ;


for i in range(len(cv2p)):
    Q_ln2k[i]=xm+cv2p[i]*dst; #caudales segun Dist. LOG PEARSON 2P
