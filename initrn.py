from numpy import *
from numpy.linalg import *

set_printoptions(precision=16)
from pyshtools.expand import SHGLQ #SHExpandDH,SHExpandDHC
#from pyshtools.rotate import djpi2,SHRotateRealCoef,SHRotateCoef
#from pyshtools.shio import SHCilmToCindex,SHCindexToCilm
from pyshtools.legendre import*

Nlr=63

Nmr=Nlr
Nch=Nlr+1

NN=Nlr*(1+(Nlr+1)//2)

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
st1=zeros((Nch,NN),dtype=float64)
mfs=zeros((Nch,NN),dtype=float64)
mft=zeros((Nch,NN),dtype=float64)
mft1=zeros((Nch,NN),dtype=float64)

for i in range(Nch):
    x,y=PlmON_d1(Nlr, xch[i])
    m0[i,:]=x[kin[:]]
    mfs[i,:]=mk[:]*x[kin[:]]/sqrt(1-xch[i]**2)
    s0[i,:]=-(nk[:]+1)*nk[:]*x[kin[:]]
    st1[i,:]=y[kin[:]]*(1-xch[i]**2)
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

m00=zeros((NN,NN),dtype=float64)
ms0=zeros((NN,NN),dtype=float64)
SS0=s0
MFT=mft
MFS=mfs

for i in range(NN):
    ji=where(mk[i] == mk[range(NN)])[0]
    for j in range(NN): # ji:
        if mk[i] == mk[j]:
            ms0[i,j]=sum(m0[:,i]*s0[:,j]*wch[:])
ms_0=inv(ms0)

    
    
MRT=zeros((NN,NN),dtype=float64)
mss=zeros((NN,NN),dtype=float64)
ms1=zeros((NN,NN),dtype=float64)
m0m=zeros((NN,NN),dtype=float64)
mpp=zeros((NN,NN),dtype=float64)
m1ss=zeros((NN,NN),dtype=float64)


for i in range(NN):
    ji=where(mk[i] == mk[range(NN)])[0]
    for j in range(NN): # ji:
        if mk[i] == mk[j]:
            m00[i,j]=sum(m0[:,i]*m0[:,j]*wch[:])
            MRT[i,j]=sum(m0[:,i]*mk[j]*m0[:,j]*(0.-xch[:]**2/4.)*wch[:])
            mss[i,j]=sum(m0[:,i]*mk[j]*s0[:,j]*(0.-xch[:]**2/4.)*wch[:])
            ms1[i,j]=sum(m0[:,i]*st1[:,j]*wch[:])
            m0m[i,j]=sum(m0[:,i]*mk[j]**2*m0[:,j]*wch[:])#/(1-xch**2))
            mpp[i,j]=sum(m0[:,i]*(s0[:,j]*xch[:]**2)*wch[:])
            m1ss[i,j]=sum(m0[:,i]*mk[j]**2*s0[:,j]*wch[:])

MB1=ms0 #+4*msm #dot(ms_0,ms1)
MRS=mss
MSA=dot(ms_0,m1ss)
MM0=m0
M0M=m0m
#
for i in range(NN):
    for j in range(NN): # ji:
        if mk[i] == mk[j]:
            if mk[i]==0:
                fac=2*pi
            else:
                fac=pi
            m00[i,j]=sum(m0[:,i]*m0[:,j]*wch[:])*fac
            ms0[i,j]=sum(m0[:,i]*s0[:,j]*wch[:])*fac
            msf[i,j]=sum(m0[:,i]*mfs[:,j]*wch[:])*fac
            mst[i,j]=sum(m0[:,i]*mft[:,j]*wch[:])*fac
            mst1[i,j]=sum(m0[:,i]*(mft[:,j]-mft1[:,j])*wch[:])*fac
#
#MM0=m0
ms_0=inv(ms0)
