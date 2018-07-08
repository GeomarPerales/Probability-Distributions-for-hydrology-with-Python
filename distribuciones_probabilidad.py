#hidrologia estadistica
import numpy as np
import math
import astropy
from astropy.table import Table
import pandas as pd
from csv import reader
import statistics as stats

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

#---------------------------------------------------------------------------------------------------
#FACTOR DE FRECUENCIA usando WEIBULL
#---------------------------------------------------------------------------------------------------
#muestra 1:59
N1=np.arange(1,60)

#coeficientes
co=2.515517
c1=0.802853
c2=0.010328

d1=1.432788
d2=0.189269
d3=0.001308

Px1=np.zeros(len(N1))
Px11=np.zeros(len(N1))
TM1=np.zeros(len(N1))
W1=np.zeros(len(N1))
k1=np.zeros(len(N1))

for i in range(len(N1)):
    Px1[i]=(i+1)/(len(N1)+1)
    Px11[i]=(i+1)/(len(N1)+1)
    TM1[i]=1/(1-Px1[i])
    if Px1[i]>=.5:
        Px1[i]=1-Px1[i]
    W1[i]=((math.log(1/Px1[i]**2)))**0.5
    k1[i]=W1[i]-(co+c1*W1[i]+c2*W1[i]**2)/(1+d1*W1[i]+d2*W1[i]**2+d3*W1[i]**3);

FF59=Table([TM1,N1,Px1,W1,k1])



#muestra 1-99
N2=np.arange(1,100)

Px2=np.zeros(len(N2))
Px22=np.zeros(len(N2))
TM2=np.zeros(len(N2))
W2=np.zeros(len(N2))
k2=np.zeros(len(N2))

for i in range(len(N2)):
    Px2[i]=(i+1)/(len(N2)+1)
    Px22[i]=(i+1)/(len(N2)+1)
    TM2[i]=1/(1-Px2[i])
    if Px2[i]>=.5:
        Px2[i]=1-Px2[i]
    W2[i]=((math.log(1/Px2[i]**2)))**0.5
    k2[i]=W2[i]-(co+c1*W2[i]+c2*W2[i]**2)/(1+d1*W2[i]+d2*W2[i]**2+d3*W2[i]**3);

FF99=Table([TM2,N2,Px2,W2,k2])



#muestra 1-199
N3=np.arange(1,200)

Px3=np.zeros(len(N3))
Px33=np.zeros(len(N3))
TM3=np.zeros(len(N3))
W3=np.zeros(len(N3))
k3=np.zeros(len(N3))

for i in range(len(N3)):
    Px3[i]=(i+1)/(len(N3)+1)
    Px33[i]=(i+1)/(len(N3)+1)
    TM3[i]=1/(1-Px3[i])
    if Px3[i]>=.5:
        Px3[i]=1-Px3[i]
    W3[i]=((math.log(1/Px3[i]**2)))**0.5
    k3[i]=W3[i]-(co+c1*W3[i]+c2*W3[i]**2)/(1+d1*W3[i]+d2*W3[i]**2+d3*W3[i]**3);

FF199=Table([TM3,N3,Px3,W3,k3])


            
#muestra 1-499
N4=np.arange(1,500)

Px4=np.zeros(len(N4))
Px44=np.zeros(len(N4))
TM4=np.zeros(len(N4))
W4=np.zeros(len(N4))
k4=np.zeros(len(N4))

for i in range(len(N4)):
    Px4[i]=(i+1)/(len(N4)+1)
    Px44[i]=(i+1)/(len(N4)+1)
    TM4[i]=1/(1-Px4[i])
    if Px4[i]>=.5:
        Px4[i]=1-Px4[i]
    W4[i]=((math.log(1/Px4[i]**2)))**0.5
    k4[i]=W4[i]-(co+c1*W4[i]+c2*W4[i]**2)/(1+d1*W4[i]+d2*W4[i]**2+d3*W4[i]**3);

FF499=Table([TM4,N4,Px4,W4,k4])
            

       
#muestra 1-999
N5=np.arange(1,1000)

Px5=np.zeros(len(N5))
Px55=np.zeros(len(N5))
TM5=np.zeros(len(N5))
W5=np.zeros(len(N5))
k5=np.zeros(len(N5))

for i in range(len(N5)):
    Px5[i]=(i+1)/(len(N5)+1)
    Px55[i]=(i+1)/(len(N5)+1)
    TM5[i]=1/(1-Px5[i])
    if Px5[i]>=.5:
        Px5[i]=1-Px5[i]
    W5[i]=((math.log(1/Px5[i]**2)))**0.5
    k5[i]=W5[i]-(co+c1*W5[i]+c2*W5[i]**2)/(1+d1*W5[i]+d2*W5[i]**2+d3*W5[i]**3);

FF199=Table([TM5,N5,Px5,W5,k5])

#DISTRIBUCIONES DE PROBABILIDAD
#DISTRIBUCION NORMAL

data=pd.read_csv('data.csv')
a=data['registro']
Q=data['caudal']

xm=stats.mean(Q)
dst=stats.stdev(Q)
cv=dst/xm

TR=[2,5,10,20,50,100,200,500,1000] #tiempo de retorno
Px=[Px22[49],Px22[79],Px22[89],Px22[94],Px22[97],Px22[98],Px33[198],Px44[498],Px55[998]]
VAE=[k2[49],k2[79],k2[89],k2[94],k2[97],k2[98],k3[198],k4[498],k5[998]] #var. aleatoria standard

V_DN=[xm,dst,cv]

Q_n=np.zeros((len(VAE)))

for i in range (len(VAE)):
    Q_n[i]=xm+VAE[i]*dst #caudales segun Dist. NORMAL



        
#TABLA DE FACTOR DE FRECUENCIA para prob. 50%, 80%, 90% y 100% y
#coeficiente de variación (CV) que varian de 0.05 hasta 1.00.

cvn=np.arange(0.05,1.05,0.05)

kn=np.zeros((len(cvn),len(VAE)))

for i in range (len(cvn)):
    for j in range(len(VAE)):
        kn[i][j]=(math.exp((math.log(1+cvn[i]**2))**0.5*VAE[j]-0.5*(math.log(1+cvn[i]**2)))-1)/cvn[i]



#LOG NORMAL 2P
lnQ=np.log(Q)
cvQ=dst/xm

uy1=0.5*math.log(xm**2/(1+cvQ**2))
dy1=(math.log(1+cvQ**2))**0.5

V_DN2P=[cvQ,uy1,dy1]

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


#LOG NORMAL 3P
sg=0

#coeficiente asimetria
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



#TABLA DE FACTOR DE FRECUENCIA PARA GUMBEL
#PARA PROB. DE 50%,80%,90%, 95% Y 100 CON INC. DE 5 U

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

for i in range(len(TR)):
    f_TR[i]=-np.log(-np.log((TR[i]-1)/TR[i]))

Yt=np.zeros((len(TM),len(f_TR)))

for i in range(len(TM)):
    for j in range(len(f_TR)):
        Yt[i][j]=-(M[i]-f_TR[j])/DS[i]
    


#DISTRIBUCION GUMBEL
al=1.2825/dst #alfa
u=xm-0.45*dst #u

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


#TABLA DE FRECUENCIA PEARSON
#prob .... y coe;f de sesgo de 0.0 a 2.0 con inc de 0.1
csp=np.arange(0,2.1,0.1)

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


#DISTRIBUCION PEARSON
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


#DISTRIBUCIÓN LOG PEARSON
xm_lp=stats.mean(lnQ)
ds_lp=stats.stdev(lnQ)
cv_lp=ds_lp/xm_lp

V_LP=[xm,ds_lp,cv_lp]

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

