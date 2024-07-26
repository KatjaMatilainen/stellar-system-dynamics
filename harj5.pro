common h5_common,G,M,a,v0,rc,q

;Vakiot ja alkuarvot
x0=0.6d0
y0=0.d0
vx0=0.5d0
v0=1.d0
rc=0.14d0
q=0.9d0
E0=-0.337d0
pot0=0.5d0*v0^2*alog(rc^2+x0^2+y0^2/q^2)
vy0=sqrt((E0-pot0)*2-vx0^2)
;vy0=0.4d0
G=1.d0
M=1.d0
a=1.d0
T=2.*!dpi
dt=T/1000
nstep=10*T/dt

;Vektorit lopputuloksille
x_vec=findgen(nstep)*0.d0
x_vec[0]=x0
y_vec=x_vec
y_vec[0]=y0
vx_vec=x_vec
vx_vec[0]=vx0
vy_vec=x_vec
vy_vec[0]=vy0
timet=findgen(nstep)*dt

;RK4-looppi
time=0
yy=[x0,y0,vx0,vy0]
for i=0,nstep-1 do begin
  der=deriv_nsym(time,yy)
  res=rk4(yy,der,time,dt,'deriv_nsym',/double)
  yy=res

  x_vec[i]=yy[0]
  y_vec[i]=yy[1]
  vx_vec[i]=yy[2]
  vy_vec[i]=yy[3]
endfor

;pot=0.5d0*v0^2*alog(rc^2+x_vec^2+y_vec^2/q^2)
;E=pot+0.5d0*(vx_vec^2+vy_vec^2)
;L=x_vec*vy_vec-y_vec*vx_vec

;!p.multi=[0,2,2]
nwin
;plot,x_vec,y_vec,/iso
;print,x_vec,y_vec
;plot,timet,E-E[0],title='Energia',xr=[0,60]
;plot,timet,L

;Contour
xc=-2.+findgen(401)*0.01d0
yc=-2.+findgen(401)*0.01d0
n=n_elements(xc)
pot=fltarr(401,401)

for i=0,n-1 do begin
  for j=0,n-1 do begin
     pot(i,j)=0.5d0*v0^2*alog(rc^2+xc[i]^2+yc[j]^2/q^2)
  endfor
endfor

contour,pot,xc,yc,nlev=10,c_col=lindgen(10)+2,/iso
contour,pot,xc,yc,lev=E0,/over,thick=2,col=2
oplot,x_vec,y_vec
;print,pot

x_sos=x_vec*0.d0
vx_sos=x_sos
n_sos=0
m=n_elements(y_vec)

;Surface of section
for i=1,m-1 do begin
  if (y_vec[i-1] lt 0 and y_vec[i] gt 0) then begin
     x_sos[n_sos]=x_vec[i] 
     vx_sos[n_sos]=vx_vec[i]
     n_sos=n_sos+1
  endif
endfor

x_sos=x_sos[0:n_sos-1]
vx_sos=vx_sos[0:n_sos-1]

nwin
plot,x_sos,vx_sos,psym=6

end
