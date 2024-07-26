common h7_common,G,M,a,eps,omeb,lambda

;H7: Resonances in non-axisymmetric rotating potential
;Toomre disk

;------------------------------------------------------;
;Vakiot ja alkuarvot
;------------------------------------------------------;

;Vakiot
G=1.d0
M=1.d0
a=0.25d0
eps=0.1d0
omeb=0.75d0
T=2.d0*!dpi
rmax=2.d0
lambda=0.5d0

;Alkuarvot
n0=10000l
x0=(randomu(seed,n0,/double)*2.-1.)*rmax
y0=(randomu(seed,n0,/double)*2.-1.)*rmax
r0=sqrt(x0^2+y0^2)
f0=atan(y0,x0)
f0=f0(where (r0 le rmax))
r0=r0(where (r0 le rmax))
n=n_elements(r0)

rp0=r0*0.d0
fp0=sqrt(1.d0/(r0^2+a^2)^1.5d0)-omeb

!except=0
;stop

;Aika-askel
dt=T/200
nstep=10*T/dt

;Jacobi constant
apu_b=-eps*a*r0^2/(r0^2+a^2)^2
pot0=-1.d0/sqrt(r0^2+a^2)
pot1=apu_b*cos(2.d0*f0)

pot=pot0+pot1

ej0=pot-0.5d0*r0^2*omeb^2+0.5d0*(rp0^2+(r0*fp0)^2)

;Resonanssirajat:

;CR
x=findgen(10000)/1000
omega=(x^2+a^2)^(-0.75)
apu1=4-3*x^2*(x^2+a^2)^(-1)
apu2=(x^2+a^2)^1.5
kappa=sqrt(apu1/apu2)

test_cr=omega-omeb
ind_cr=where(test_cr lt 0)
rcr=x(ind_cr(0))
print,'r_cr=',rcr

;ILR
test_ilr=omega-kappa/2-omeb
ind_ilr=where((test_ilr lt 0) and (x gt 0.5))
rilr=x(ind_ilr(0))
print,'r_ilr',rilr

;OLR
test_olr=omega+kappa/2-omeb
ind_olr=where(test_olr lt 0)
rolr=x(ind_olr(0))
print,'r_olr',rolr

nwin
plot,x,omega
oplot,x,omega-kappa/2
oplot,x,omega+kappa/2
oplot,x,x*0+omeb

;--------------------------------------------------------;
; RK4
;--------------------------------------------------------;


nwin,xs=400,ys=400

time=0.
yy=[r0,f0,rp0,fp0]
for i=0,nstep-1 do begin
  der=deriv_h7(time,yy)
  res=rk4(yy,der,time,dt,'deriv_h7',/double)
  yy=res

r=yy[0l:n-1]
f=yy[n:2*n-1]
rp=yy[2*n:3*n-1]
fp=yy[3*n:4*n-1]


plot,r*cos(f),r*sin(f),psym=3,xr=[-1,1]*2.*rmax,yr=[-1,1]*2.*rmax,title=lambda
oplot,rcr*cos(f),rcr*sin(f),psym=3,col=3
oplot,rilr*cos(f),rilr*sin(f),psym=3,col=4
oplot,rolr*cos(f),rolr*sin(f),psym=3,col=5

;Bitmap:
;loadct,0,/sil
;xy_to_dens3,x,y,400,400,rmax*2,rmax*2,dens,/nop,/nor,sil=3
;tvscl,alog(dens>.1)

endfor

;Jacobi constant
apu_b=-eps*a*r^2/(r^2+a^2)^2
pot0=-1.d0/sqrt(r^2+a^2)
pot1=apu_b*cos(2.d0*f)

pot=pot0+pot1

ej=pot-0.5d0*r^2*omeb^2+0.5d0*(rp^2+(r*fp)^2)

nwin
plot,r0,(ej-ej0),psym=3,yr=[-0.05,0.05]*0

end
