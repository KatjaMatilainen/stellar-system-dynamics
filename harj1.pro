;P‰‰ohjelma, joka plottaa erilaisten pallosymmetristen tiheysjakaumien 
;massan,voiman,ratanopeuden jne. s‰teen funktiona.

;Yleiset asetukset:
ps=0
a=1.d0
G=1.
b=1.
fac=1.
xr=[0,5]
yr=[0,10]
imagemax=1./50./50.*10

;--------------------------------------------------------------------;
; B) Plummer sphere
;--------------------------------------------------------------------;

r=1.d-6*1.01d0^dindgen(2000)
index=n_elements(r)
rho=findgen(index)*0.d0

;-------------------------------------------------------;
;Kutsutaan aliohjelmaa plummer.pro
plummer,r,rho=rho,mass=mass,grav=grav,pot=pot,vrot=vrot, $
vesc=vesc,omega=omega,kappa=kappa
;-------------------------------------------------------;

;plotataan malli

psdirect,'harj1_plummer',ps,/color
!p.multi=[0,2,2]
nwin

;(X,Y)=(r,rho)
plot,r,rho,xtitle='r',ytitle='rho(r)',xr=xr,psym=-6,syms=.1,yr=yr,$
title='Plummer sphere'

;(X,Y)=(r,vrot(r))
plot,r,vrot,xtitle='r',ytitle='vrot(r)',xr=xr,yr=[0,5],psym=-6,syms=.1

;---------------------------------------------------------------------;
;Lasketaan teoreettiset arvot potentiaalille ja ratanopeudelle
;---------------------------------------------------------------------;
M=max(mass)

pot_teor=-G*M*(r^2+b^2)^(-0.5)

vcirc_teor=sqrt(G*M*r^2/(r^2+b^2)^1.5)

;--------------------------------------------------------------------;
;Plotataan vertailun vuoksi samaan kuvaan teor. ja numeeriset suureet
;--------------------------------------------------------------------;

plot,r,pot_teor,xtitle='r',ytitle='pot(r)', xr=xr,yr=[-30,1],col=2
oplot,r,pot,psym=-6,syms=.1

plot,r,vcirc_teor,xr=xr,yr=yr,xtitle='r',ytitle='v_rot(r)'
oplot,r,vrot,col=2

;------------------------------------------------;
; HOMOGENOUS SPHERE
;------------------------------------------------;

;Teoreettiset arvot:
ind1=where(r lt a)
ind2=where(r gt a)

Mtot=4./3*!dpi*a^3

vrot_h_teor=r*0.d0
vrot_h_teor(ind1)=sqrt(G*Mtot/a)*r(ind1)/a
vrot_h_teor(ind2)=sqrt(G*Mtot/r(ind2))

;-----------------------------------------------------------------;
;Kutsutaan ohjelmaa homogenous.pro
homogenous,r,rho=rho_h,mass=mass_h,grav=grav_h,pot=pot_h,vrot=vrot_h
;-----------------------------------------------------------------;

;Plotataan teor. ja numeer. arvot samaan kuvaan
!p.multi=[0,2,2]
nwin
plot,r,vrot_h_teor,xtitle='r',ytitle='v_rot(r)',xr=xr,yr=yr,$
title='Homogenous sphere'
oplot,r,vrot_h,col=2

;Potentiaali
plot,r,pot_h,xr=xr,yr=[-10,1],xtitle='r',ytitle='pot(r)'

;Tiheysjakauma
plot,r,rho_h,xr=[0,2],yr=[0,2],xtitle='r',ytitle='rho(r)'

;-------------------------------------------------------------------;
; ISOTERMINEN
;-------------------------------------------------------------------;
rc=1.d0
rh=10.d0
Mtot_i=4.*!dpi*rh*rc^2

;Teoreettiset arvot
pot_i_teor=G*Mtot_i/(2*rh)*alog(r^2+rc^2)
pot_i_teor=pot_i_teor-max(pot_i_teor)

vrot_i_teor=sqrt(G*Mtot_i/rh*r^2/(r^2+rc^2))

print,'vrot_i_teor',vrot_i_teor

;-------------------------------------------------------------------;
;Kutsutaan ohjelmaa isothermal.pro
isothermal,r,rho=rho_i,mass=mass_i,grav=grav_i,pot=pot_i,vrot=vrot_i
;-------------------------------------------------------------------;

;Plotataan teor. ja num. arvot

!p.multi=[0,2,2]
nwin

plot,r,pot_i_teor,xr=[0,30],yr=[-100,1],xtitle='r',ytitle='pot(r)',$
title='Isothermal sphere'
oplot,r,pot_i,col=2

plot,r,vrot_i_teor,xr=xr,yr=yr,xtitle='r',ytitle='vrot(r)'
oplot,r,vrot_i,col=2

end
