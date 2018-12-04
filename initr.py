from numpy import *
from numpy.linalg import *
from scipy.special import *

set_printoptions(precision=16)
from pyshtools.expand import SHGLQ #SHExpandDH,SHExpandDHC
#from pyshtools.rotate import djpi2,SHRotateRealCoef,SHRotateCoef
#from pyshtools.shio import SHCilmToCindex,SHCindexToCilm
from pyshtools.legendre import*

Nlr=299

Nmr=Nlr
Nch=Nlr+1

k=0
for j in range(0,Nmr+1):
    if j==0:
        for i in range(j+1,Nlr+1): #range(-i-1,i+2): # range(1): #range(-1,2):
            k=k+1
    else:
        for i in range(j,Nlr+1): #range(-i-1,i+2): # range(1): #range(-1,2):
            k=k+1


NN=k #Nlr*Nmr #int((Nn+1m)*Nn/2)

k=0
for j in range(0,Nmr+1):
    if j==0:
        for i in range(j+1,Nlr+1): #range(-i-1,i+2): # range(1): #range(-1,2):
            k=k+1
    else:
        for i in range(j,Nlr+1): #range(-i-1,i+2): # range(1): #range(-1,2):
            k=k+1
NN0=k #Nlr*Nmr #int((Nn+1m)*Nn/2)

xchl,wchl=SHGLQ(Nlr)


xch=((xchl[:]-xchl[-1::-1])/2.)[-1::-1]
wch=((wchl[:]+wchl[-1::-1])/2.)[-1::-1]

xchl1,wchl1=SHGLQ(Nlr-1)
xch1=((xchl1[:]-xchl1[-1::-1])/2.)[-1::-1]



mk=zeros(NN,dtype=int)
mkk=zeros(NN,dtype=int)
nk=zeros(NN,dtype=int)
kin=zeros(NN,dtype=int)

k=0
for j in range(0,Nmr+1):
    if j==0:
        for i in range(j+1,Nlr+1): 
            nk[k]=i  
            mk[k]=j
            kin[k]=nk[k]*(nk[k]+1)//2+mk[k]
            k=k+1
    else:
        for i in range(j,Nlr+1): 
            nk[k]=i  
            mk[k]=j
            kin[k]=nk[k]*(nk[k]+1)//2+mk[k]
            k=k+1


m0=zeros((Nch,NN),dtype=float64)
m_0=zeros((Nch,NN),dtype=float64)
m1=zeros((Nch,NN),dtype=float64)
s0=zeros((Nch,NN),dtype=float64)
mfs=zeros((Nch,NN),dtype=float64)
mft=zeros((Nch,NN),dtype=float64)
mft1=zeros((Nch,NN),dtype=float64)

for i in range(Nch):
    x,y=PlmON_d1(Nlr, xch[i])
    m0[i,:]=x[kin[:]]
    mfs[i,:]=mk[:]*x[kin[:]]/sqrt(1-xch[i]**2)
    s0[i,:]=-(nk[:]+1)*nk[:]*x[kin[:]]
    mft[i,:]=y[kin[:]]*sqrt(1-xch[i]**2)
    mft1[i,:]=x[kin[:]]*xch[i]/sqrt(1-xch[i]**2)

m00=zeros((NN,NN),dtype=float64)
ms0=zeros((NN,NN),dtype=float64)
SS0=s0
MFT=mft
MFS=mfs
msf=zeros((NN,NN),dtype=float64)
mst=zeros((NN,NN),dtype=float64)
mst1=zeros((NN,NN),dtype=float64)

for i in range(NN):
    for j in range(NN): # ji:
        if mk[i] == mk[j]:
            m00[i,j]=sum(m0[:,i]*m0[:,j]*wch[:])*2*pi
            ms0[i,j]=sum(m0[:,i]*s0[:,j]*wch[:])*2*pi
            msf[i,j]=sum(m0[:,i]*mfs[:,j]*wch[:])*2*pi
            mst[i,j]=sum(m0[:,i]*mft[:,j]*wch[:])*2*pi
            mst1[i,j]=sum(m0[:,i]*(mft[:,j]-mft1[:,j])*wch[:])*2*pi

MM0=m0
ms_0=inv(ms0)
