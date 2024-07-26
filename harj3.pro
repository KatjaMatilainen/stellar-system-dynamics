;---------------------------------------------------------------------;
; Harj.3) Ratojen numeerinen integrointi
;---------------------------------------------------------------------;

;3.1) Toomre-Kuzmin-levyn ekvaattoritaso
; I asteen Taylor

;--------------------------------------------------; 
;Vakiot
x0=1.d0
y0=0.d0
vx0=0.d0
vy0=1.d0
G=1.d0
M=1.d0
a=0.1d0
ax0=-G*M*x0*(x0^2+y0^2+a^2)^(-3./2)
ay0=-G*M*y0*(x0^2+y0^2+a^2)^(-3./2)

T=2.*!dpi
dt=T/1000
nstep=10*T/dt
;-------------------------------------------------;

x_vec=findgen(nstep)*0.d0
x_vec[0]=x0
vx_vec=x_vec
vx_vec[0]=vx0
y_vec=x_vec
y_vec[0]=y0
vy_vec=x_vec
vy_vec[0]=y0
ax_vec=x_vec
ax_vec[0]=ax0
ay_vec=x_vec
ay_vec[0]=ay0
timet=findgen(nstep)*dt

for i=1,nstep-1 do begin
  ax0=-G*M*x0*(x0^2+y0^2+a^2)^(-3./2)
  ay0=-G*M*y0*(x0^2+y0^2+a^2)^(-3./2)
  x0=x0+vx0*dt
  y0=y0+vy0*dt
  vx0=vx0+ax0*dt
  vy0=vy0+ay0*dt

;Tallennetaan vektoreihin
  x_vec[i]=x0
  y_vec[i]=y0
  vx_vec[i]=vx0
  vy_vec[i]=vy0
  ax_vec[i]=ax0
  ay_vec[i]=ay0
endfor

pot=-G*M/(sqrt(x_vec^2+y_vec^2+a^2))
E=pot+0.5d0*(vx_vec^2+vy_vec^2)
L=x_vec*vy_vec-y_vec*vx_vec

!p.multi=[0,2,2]
nwin
plot,x_vec,y_vec,title=a
plot,timet,E,title='Energia  (I Taylor)',xr=[0,60],yr=[-0.55,-0.3]
plot,timet,L,title='Kulmaliikem‰‰r‰  (I Taylor)',xr=[0,60],yr=[0.98,1.23]

end
