common h6_common,G,M,a,v0,rc,re,q

;Vakiot ja alkuarvot
vx0=0.1d0
vy0=0.d0
v0=1.d0
rc=0.1d0
re=1.d0
q=0.9d0
E0=-0.337d0
x0=0.8
y0=0.d0
r0=sqrt(x0^2+y0^2)
phi0=atan(y0/x0)
pot0=0.5d0*v0^2*alog(rc^2+x0^2+y0^2/q^2-r0^3/re*cos(2*phi0))
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

;RK4
time=0
yy=[x0,y0,vx0,vy0]
for i=0,nstep-1 do begin
   der=deriv_rot(time,yy)
   res=rk4(yy,der,time,dt,'deriv_rot',/double)
   yy=res
   
  x_vec[i]=yy[0]
  y_vec[i]=yy[1]
  vx_vec[i]=yy[2]
  vy_vec[i]=yy[3]
endfor

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

end
