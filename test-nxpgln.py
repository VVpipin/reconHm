from numpy import *
from numpy.fft import *
from scipy import *
from scipy.interpolate import *
from pylab import *


from scipy.integrate import *
    
import matplotlib.cm as cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

params = {'axes.labelsize': 24,
          'font.size': 24,
          'legend.fontsize': 30,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'savefig.dpi':32,
          'text.usetex': True}
plt.rcParams.update(params)
from pyshtools import *
from pyshtools.expand import SHGLQ,SHExpandGLQ,MakeGridGLQ,SHExpandDH,MakeGridDH,SHMultiply
#from pyshtools.rotate import djpi2,SHRotateRealCoef,SHRotateCoef
from pyshtools.shio import SHCilmToCindex,SHCindexToCilm,SHrtoc,SHctor
from pyshtools.spectralanalysis import spectrum,cross_spectrum,SHAdmitCorr,SHConfidence 


npzfile=load('bm.npz')
br=npzfile['br']
bf=npzfile['bf']
bp=npzfile['bt']

npzfile=load('st.npz')
ss=npzfile['s']
tt=npzfile['t']

Nph0=72
phi0=arange(Nph0)*2*pi/Nph0

from initrn import *
Nphg=Nlr*2+1
phig=arange(Nphg)*2*pi/Nphg

f = RectBivariateSpline(xch,phi0,br)
brn = f(xch,phig)

f = RectBivariateSpline(xch,phi0,bf)
bfn = f(xch,phig)

f = RectBivariateSpline(xch,phi0,bp)
bpn = f(xch,phig)
zero=xch[-1::-1]

xg,yg=meshgrid(phig*180/pi,90-arccos(zero)*180/pi)
cnilm=zeros((2,Nlr+1,Nmr+1))


#################
#################
brr=brn#tensordot(mean(brn,axis=1),ones(Nphg),0)
cilm = SHExpandGLQ(brr[-1::-1,:],wch,zero, norm=4)
ccilm=SHrtoc(cilm)
brco=(ccilm[0,nk,mk]+1j*ccilm[1,nk,mk]);
s0r=-dot(ms_0,brco)
ath=1j*dot(msf,s0r);aff=dot(mst,s0r)    

cilm = SHExpandGLQ(bfn[-1::-1,:],wch,zero, norm=4)
ccilm=SHrtoc(cilm)
bfco=ccilm[0,nk,mk]+1j*ccilm[1,nk,mk];

cilm = SHExpandGLQ(bpn[-1::-1,:],wch,zero, norm=4)
ccilm=SHrtoc(cilm)
bpco=(ccilm[0,nk,mk]+1j*ccilm[1,nk,mk]);

t0r=dot(ms_0, dot(mst1,bfco)+1j*dot(msf,bpco))
fsr=dot(ms_0, 1j*dot(msf,bfco)-dot(mst1,bpco))

cilm=zeros((2,Nlr+1,Nmr+1))
brcc=-dot(ms0,s0r)
cilm[0,nk,mk]=real(brcc);cilm[1,nk,mk]=imag(brcc) 
rcilm=SHctor(cilm)
brc=(rcilm)
brnr = MakeGridGLQ(rcilm,zero,norm=4)

cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(s0r);
cilm[1,nk,mk]=imag(s0r) 
src=SHctor(cilm)
srm = MakeGridGLQ(src,zero,norm=4) #,lmax_calc=Nlr-10)

cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(t0r);
cilm[1,nk,mk]=imag(t0r) 
arc=SHctor(cilm)
ar = MakeGridGLQ(arc,zero,norm=4) 


cilm = SHExpandGLQ(ar,wch,zero, norm=4)
ccilm=SHrtoc(cilm)
curb=-dot(ms0,ccilm[0,nk,mk]+1j*ccilm[1,nk,mk]);
cilm[0,nk,mk]=real(curb);
cilm[1,nk,mk]=imag(curb); 
curlc=SHctor(cilm)
curlr = MakeGridGLQ(curlc,zero,norm=4)


cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(fsr);
cilm[1,nk,mk]=imag(fsr) 
fsrc=SHctor(cilm)
fsrm = MakeGridGLQ(fsrc,zero,norm=4) #,lmax_calc=Nlr-10)

   
cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(ath);cilm[1,nk,mk]=imag(ath) 
rcilm=SHctor(cilm)
apc=(rcilm)
ap = MakeGridGLQ(apc,zero,norm=4)#,lmax_calc=Nlr*2//3)

cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(bpco);cilm[1,nk,mk]=imag(bpco) 
bpc=SHctor(cilm)
bpnr = MakeGridGLQ(bpc,zero,norm=4)#,lmax_calc=Nlr*2//3)


cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(aff);cilm[1,nk,mk]=imag(aff) 
afc=SHctor(cilm)
af = MakeGridGLQ(afc,zero,norm=4)#,lmax_calc=Nlr*2//3)


cilm=zeros((2,Nlr+1,Nmr+1))
cilm[0,nk,mk]=real(bfco);cilm[1,nk,mk]=imag(bfco) 
bfc=SHctor(cilm) 
bfnr = MakeGridGLQ(bfc,zero,norm=4)#,lmax_calc=Nlr*2//3)


   
hln=ar*brnr+af*bfnr+ap*bpnr
hc=brnr*curlr


npzfile=load('hh.npz')
hcm=npzfile['hc']
hmm=npzfile['hm']

f = RectBivariateSpline(xch,phi0,hcm)
hcn = f(xch,phig)[-1::-1,:]
f = RectBivariateSpline(xch,phi0,hmm)
hmn = f(xch,phig)[-1::-1,:]

#there is a small systematic difference between results which is 
# because of interpolation and different decomposition schemes employed in
#SHTools and in my dynamo model where I work with truncated series over m
# using the "romb" scheme (m_max=l_max/2) while SHTools use the more accurate decomposition
# where m_max==l_max  
#below are figures for Br,AB,hc and S all of them in well agreement with 
#originals br,hmn,hcn and ss

plt.figure()
CS2=plt.pcolor(xg,yg,brnr,cmap=cmap.seismic,vmin=-.2,vmax=.2)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB=plt.colorbar(CS2, cax=cax,ticks=[-.1,0,.1])
CB.set_label('B$_{r}$',fontsize=30)
plt.tight_layout()
plt.show()

plt.figure()
CS2=plt.pcolor(xg,yg,hln*10,cmap=cmap.seismic,vmin=-.2,vmax=.2)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB=plt.colorbar(CS2, cax=cax,ticks=[-.1,0,.1])
CB.set_label('10AB',fontsize=30)
plt.tight_layout()
plt.show()


plt.figure()
CS2=plt.pcolor(xg,yg,100*srm,cmap=cmap.seismic,vmin=-.2,vmax=.2)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB=plt.colorbar(CS2, cax=cax,ticks=[-.1,0,.1])
CB.set_label('100S',fontsize=30)
plt.tight_layout()
plt.show()

plt.figure()
CS2=plt.pcolor(xg,yg,hc,cmap=cmap.seismic,vmin=-.2,vmax=.2)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB=plt.colorbar(CS2, cax=cax,ticks=[-.1,0,.1])
CB.set_label('hc',fontsize=30)
plt.tight_layout()
plt.show()
