#import sys
fim='hmi.B_synoptic_small.'
fima='/home/va/work/HMI/'
#sys.path.append(fim);
fnam='hmi.B_synopticI.fits'
from numpy import *
from numpy.fft import *
from scipy import *
from scipy.interpolate import *
from pylab import *


from astropy.io import fits
from scipy.integrate import *
    
import matplotlib.cm as cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.utils.data import get_pkg_data_filename
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


CR=range(2097,2181)
LM=len(CR)
NT0=360
Nph0=720

Rs=6.96e10
mu0=linspace(-1,1,NT0)
phi0=arange(Nph0)*2*pi/Nph0


c1=sum(ones(Nph0)*simps(ones(NT0),mu0))*2*pi/Nph0
c1=simps(ones(Nph0)*simps(ones(NT0),mu0),phi0)
LM=len(CR)


from initr import *
Nphg=Nlr*2+1
phig=arange(Nphg)*2*pi/Nphg
NT=2*(Nch+1);Nph=2*NT #
mu=cos(arange(NT0)*pi/(NT0*1.))[-1::-1];phi=arange(Nph0)*2*pi/Nph0

zero=xch[-1::-1]

x0,y0=meshgrid(phi0*180/pi,90-arccos(mu0)*180/pi)
xg,yg=meshgrid(phig*180/pi,90-arccos(zero)*180/pi)
brM=zeros((Nch,len(CR)))
bfM=zeros((Nch,len(CR)))
hM=zeros((Nch,len(CR)))
hMr=zeros((Nch,len(CR)))
afM=zeros((Nch,len(CR)))
arM=zeros((Nch,len(CR)))
sM=zeros((Nch,len(CR)))
hLt=zeros((Nch,len(CR)))
hS=zeros((Nch,len(CR)))
hR=zeros((Nch,len(CR)))
hLT=zeros((Nch,len(CR)))
hC=zeros((Nch,len(CR)))
hCr=zeros((Nch,len(CR)))
hcL=zeros((Nch,len(CR)))
hsp=zeros((Nlr+1,len(CR)))
hcsp=zeros((Nlr+1,len(CR)))
hhc=zeros((Nlr+1,len(CR)))
hhcr=zeros((Nlr+1,len(CR)))
hlnr=zeros(shape(xg))
hcr=zeros(shape(xg))

phin=arange(Nphg+1)*2*pi/Nphg
xn,yn=meshgrid(phin*180/pi,90-arccos(zero)*180/pi)
zero1,wch1=SHGLQ(2*Nlr)
Nlr1=2*Nlr
Nphg1=Nlr1*2+1
phig1=arange(Nphg1)*2*pi/Nphg1
xg1,yg1=meshgrid(phig1*180/pi,90-arccos(zero1)*180/pi)

import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal

i0=0
cnilm=zeros((2,Nlr+1,Nmr+1))

for i in CR:
    imff=nan_to_num(fits.getdata(fim+str(CR[i0])+'.Br.fits',ext=0))#,kernel_size=[11,21]
    BR=imff
    f = RectBivariateSpline(mu0,phi0,imff)
    imsm = f(mu,phi)
    f = RectBivariateSpline(mu,phi,imsm)
    cilm = SHExpandDH(imsm[-1::-1,:],norm=4,sampling=2,lmax_calc=Nlr)
    cilmc = SHCoeffs.from_array(cilm)
    cilmc.set_coeffs(0,0,0)
    cilmc.set_coeffs(0,arange(Nlr*3//4,Nlr),0)
    cilmo=cilmc.convert(normalization='ortho')
    cnilm = cilmo.to_array(normalization='ortho') #, csphase, lmax]
    Z1 =MakeGridGLQ(cnilm,zero,norm=4)/sqrt(4*pi)
    imf=sum(sum(Z1,axis=1)*wch)*2*pi/Nphg/4./pi

    brn =Z1-imf
    cilm = SHExpandGLQ(brn,wch,zero, norm=4)
    ccilm=SHrtoc(cilm)
    brco=(ccilm[0,nk,mk]+1j*ccilm[1,nk,mk]);
    s0r=-dot(ms_0,brco)
    ath=1j*dot(msf,s0r);aph=dot(mst,s0r)


    imff=nan_to_num(fits.getdata(fim+str(CR[i0])+'.Bt.fits',ext=0))#,kernel_size=[11,21]

    f = RectBivariateSpline(mu0,phi0,imff)
    imsm =f(mu,phi)   
    cilm = SHExpandDH(imsm[-1::-1,:],norm=4,sampling=2,lmax_calc=Nlr)
    cilmc = SHCoeffs.from_array(cilm)
    
    cilmc.set_coeffs(0,0,0)
    cilmc.set_coeffs(0,arange(Nlr*3//4,Nlr),0)
    cilmo=cilmc.convert(normalization='ortho')
    cnilm = cilmo.to_array(normalization='ortho') #, csphase, lmax]
    Z1 =MakeGridGLQ(cnilm,zero,norm=4)/sqrt(4*pi)

    imf=sum(sum(Z1,axis=1)*wch)*2*pi/Nphg/4./pi
    bfn =Z1-imf
    cilm = SHExpandGLQ(bfn,wch,zero, norm=4)
    ccilm=SHrtoc(cilm)
    bfco=ccilm[0,nk,mk]+1j*ccilm[1,nk,mk];




    imff=nan_to_num(fits.getdata(fim+str(CR[i0])+'.Bp.fits',ext=0))#,kernel_size=[11,21]
    f = RectBivariateSpline(mu0,phi0, imff)
    imsm =f(mu,phi)
    cilm = SHExpandDH(imsm[-1::-1,:],norm=4,sampling=2,lmax_calc=Nlr)
    cilmc = SHCoeffs.from_array(cilm)
    
    cilmc.set_coeffs(0,0,0)
    cilmc.set_coeffs(0,arange(Nlr*3//4,Nlr),0)
    cilmo=cilmc.convert(normalization='ortho')
    cnilm = cilmo.to_array(normalization='ortho') #, csphase, lmax]
    Z1 =MakeGridGLQ(cnilm,zero,norm=4)/sqrt(4*pi)
    imf=sum(sum(Z1,axis=1)*wch)*2*pi/Nphg/4./pi
    bpn =Z1
    cilm = SHExpandGLQ(bpn,wch,zero, norm=4)
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
    arm = MakeGridGLQ(arc,zero,norm=4) 
    
    
    cilm = SHExpandGLQ(arm,wch,zero, norm=4)
    ccilm=SHrtoc(cilm)
    curb=-dot(ms0,ccilm[0,nk,mk]+1j*ccilm[1,nk,mk]);
    cilm[0,nk,mk]=real(curb);
    cilm[1,nk,mk]=imag(curb); 
    curlc=SHctor(cilm)
    curlr = MakeGridGLQ(curlc,zero,norm=4)/6.96e8

    
    cilm=zeros((2,Nlr+1,Nmr+1))
    cilm[0,nk,mk]=real(fsr);
    cilm[1,nk,mk]=imag(fsr) 
    fsrc=SHctor(cilm)
    fsrm = MakeGridGLQ(fsrc,zero,norm=4) #,lmax_calc=Nlr-10)

   
    cilm=zeros((2,Nlr+1,Nmr+1))
    cilm[0,nk,mk]=real(ath);cilm[1,nk,mk]=imag(ath) 
    rcilm=SHctor(cilm)
    apc=(rcilm)
    athm = MakeGridGLQ(apc,zero,norm=4)#,lmax_calc=Nlr*2//3)

    cilm=zeros((2,Nlr+1,Nmr+1))
    cilm[0,nk,mk]=real(bpco);cilm[1,nk,mk]=imag(bpco) 
    bpc=SHctor(cilm)
    bpnr = MakeGridGLQ(bpc,zero,norm=4)#,lmax_calc=Nlr*2//3)


    cilm=zeros((2,Nlr+1,Nmr+1))
    cilm[0,nk,mk]=real(aph);cilm[1,nk,mk]=imag(aph) 
    afc=SHctor(cilm)
    aphm = MakeGridGLQ(afc,zero,norm=4)#,lmax_calc=Nlr*2//3)


    cilm=zeros((2,Nlr+1,Nmr+1))
    cilm[0,nk,mk]=real(bfco);cilm[1,nk,mk]=imag(bfco) 
    bfc=SHctor(cilm) 
    bfnr = MakeGridGLQ(bfc,zero,norm=4)#,lmax_calc=Nlr*2//3)


    #hlnr=zeros(shape(xg))
    #hcr=zeros(shape(xg))

    #argm = where( (abs(brn)==maximum_filter(abs(brn),size=(3,6))) & (abs(brn) > 50. ))
    
    hln=arm*brnr+aphm*bfnr+athm*bpnr
    #hlnr[argm[0],argm[1]]=(arm[argm[0],argm[1]]*brnr[argm[0],argm[1]]
    #+aphm[argm[0],argm[1]]*bfnr[argm[0],argm[1]]+athm[argm[0],argm[1]]*bpnr[argm[0],argm[1]])

    hc=brnr*curlr
    #hcr[argm[0],argm[1]]=brnr[argm[0],argm[1]]*curlr[argm[0],argm[1]]
    #b2=sqrt(brnr**2+bfnr**2+bpnr**2)
    ##bilm = SHExpandGLQ(b2,wch,zero, norm=4)
    #bsp = spectrum (bilm, normalization='ortho', unit='per_l')
    hilm = SHExpandGLQ(hln,wch,zero, norm=4)
    hcilm = SHExpandGLQ(hc,wch,zero, norm=4)
    hsp[:,i0] = spectrum (hilm, normalization='ortho', unit='per_l')
    hcsp[:,i0] = spectrum (hcilm, normalization='ortho', unit='per_l')
    hhc[:,i0] = SHAdmitCorr (hilm, hcilm)[2]
 
    plt.figure()
    plt.title('CR='+str(CR[i0]))   
    #CS0=plt.pcolor(x0,y0,BR,cmap=cmap.seismic,vmin=-100,vmax=100)
    CS0=plt.pcolor(xg,yg,brnr,cmap=cmap.seismic,vmin=-100,vmax=100)#,alpha=0.35)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    #CB0=plt.colorbar(CS2, cax=cax0,ticks=[-1,0,1])
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-50,0,50])
    CB0.set_label('Br,[G]',fontsize=30)
    plt.tight_layout()
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig('BHMI'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],hln[20:-20,:],norm=colors.SymLogNorm(linthresh=1, linscale=.1,
            vmin=-100.,vmax=100.),cmap=cmap.seismic)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-10,0,10],format='%.1i')
    CB0.set_label('AB/R, [G$^2$]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.76, top=0.89)
    plt.savefig('h'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],aphm[20:-20,:],cmap=cmap.seismic,vmin=-2,vmax=2)#,alpha=0.35)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-1,0,1])
    CB0.set_label('A$_{\\phi}$/R, [G]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig(fima+'af'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],arm[20:-20,:],cmap=cmap.seismic,vmin=-2,vmax=2)#,alpha=0.35)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-1,0,1])
    CB0.set_label('T, [G]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig(fima+'ar'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],srm[20:-20,:],cmap=cmap.seismic,vmin=-1,vmax=1)#,alpha=0.35)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-.5,0,.5])
    CB0.set_label('S/R, [G]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig(fima+'s'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],fsrm[20:-20,:],cmap=cmap.seismic,vmin=-5,vmax=5)#,alpha=0.35)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-3,0,3])
    CB0.set_label('$\\partial rS/\\partial r$ /R, [G]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig('fs'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')

    plt.figure()
    CS0=plt.pcolor(xg[20:-20,:],yg[20:-20,:],hc[20:-20,:]*1e5,norm=colors.SymLogNorm(linthresh=.1, linscale=1,
            vmin=-10.,vmax=10.),cmap=cmap.seismic)
    divider = make_axes_locatable(plt.gca())
    cax0 = divider.append_axes("right", "3%", pad="3%")
    CB0=plt.colorbar(CS0, cax=cax0,ticks=[-1,0,1],format='%.1i')
    CB0.set_label('h$_C$,  10$^-5$[G$^2$/M]',fontsize=30)
    subplots_adjust(hspace=0,bottom=0.1, left=0.1, right=0.8, top=0.89)
    plt.savefig('hc'+str(i0)+'.png',format='png',dpi=150);
    plt.close('all')


    brM[:,i0]=mean(brnr,axis=1)[-1::-1]
    bfM[:,i0]=mean(bfnr,axis=1)[-1::-1]
    afM[:,i0]=sum(aphm,axis=1)[-1::-1]/Nphg #
    arM[:,i0]=sum(arm,axis=1)[-1::-1]/Nphg #
    hM[:,i0]=sum(hln,axis=1)/Nphg
    hC[:,i0]=sum(hc,axis=1)/Nphg
    hMr[:,i0]=sum(hlnr,axis=1)/Nphg
    hCr[:,i0]=sum(hcr,axis=1)/Nphg
    sM[:,i0]=sum(srm,axis=1)/Nphg
    hcL[:,i0]=mean(brnr,axis=1)*mean(curlr,axis=1)*1000.
    htN=sum(sum(hlnr,axis=1)[20:Nch//2]*wch[20:Nch//2])/Nphg*2*pi
    htS=sum(sum(hlnr,axis=1)[Nch//2:-20]*wch[Nch//2:-20])/Nphg*2*pi
    
    print(i, htN,htS)
    i0=i0+1

xc,yc=meshgrid(CR,90-arccos(zero)*180/pi)
levs=-1+2*arange(10)/10
xk,yk=meshgrid(CR,arange(Nlr+1))

plt.figure(figsize=(10,5))
plt.title('croscorr AB vs hc')   
CS=plt.pcolor(xk,yk,hhc,cmap=cmap.seismic,vmin=0,vmax=1)
plt.ylabel('$\\ell$')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[0,.5,1])
plt.tight_layout()
subplots_adjust(top=0.95,bottom=0.1,left=0.1,right=0.9,hspace=0.0,wspace=0.2)
plt.show()

hS=(afM*bfM+brM*arM)[-1::-1,:]
g0 = gauss_kern(1)
g1 = gauss_kern(1)

plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,hM,cmap=cmap.seismic,vmin=-20,vmax=20)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-10,0,10])
CB0.set_label('$\\overline{AB}$/R, [G$^2$]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()

plt.figure()
CS=plt.pcolor(xc,yc,hC*1e5,cmap=cmap.seismic,vmin=-5,vmax=5)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-2,0,2])
CB0.set_label('h$_C$,  10$^-5$[G$^2$/M]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.85,hspace=0.0,wspace=0.2)
plt.show()


plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,afM[-1::-1,:],cmap=cmap.seismic,vmin=-5,vmax=5)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-3,0,3])
CB0.set_label('$\\bar{A}_{\\phi}$/R, [G]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()


plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,sM[:,:],cmap=cmap.seismic,vmin=-5,vmax=5)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-3,0,3])
CB0.set_label('$\\bar{S}$/R, [G]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()

plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,arM[-1::-1,:],cmap=cmap.seismic,vmin=-1,vmax=1)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-0.5,0,0.5])
CB0.set_label('$\\bar{T}$, [G]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()

htN=zeros(len(CR))
htS=zeros(len(CR))
hsN=zeros(len(CR))
hsS=zeros(len(CR))
hT=zeros(len(CR))
hTl=zeros(len(CR))

for i in range(len(CR)):
    htN[i]=sum(hM[20:Nch//2,i]*wch[20:Nch//2])/2
    htS[i]=sum(hM[Nch//2:-20,i]*wch[Nch//2:-20])/2
    hsN[i]=sum(hS[:Nch//2,i]*wch[:Nch//2])/2
    hsS[i]=sum(hS[Nch//2:,i]*wch[Nch//2:])/2
    hT[i]=sum(hM[:,i]*wch[:])
    hTl[i]=sum(hS[:,i]*wch[:])
    

plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,hS,cmap=cmap.seismic,vmin=-5,vmax=5)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-2,0,2])
CB0.set_label('$\\bar{A}\\bar{B}$/R, [G$^2$]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()



plt.figure(figsize=(10,5))
CS=plt.pcolor(xc,yc,brM,cmap=cmap.bwr,vmin=-20,vmax=20)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="3%")
CB0=plt.colorbar(CS, cax=cax,ticks=[-10,0,10])
CB0.set_label('$\\bar{B}_r$, [G]',fontsize=30)
plt.tight_layout()
subplots_adjust(top=0.975,bottom=0.1,left=0.1,right=0.87,hspace=0.0,wspace=0.2)
plt.show()

