;---------------------------------------------------;
; Ratojen 3D-integrointi sylinterikoordinaatistossa
;---------------------------------------------------;
common h4_common,G,M,a,v0,rc,q

r0=1.d0
z0=0.d0
phi0=0.d0
dr0=0.1d0
dz0=0.1d0
dphi0=1.d0
G=1.d0
M=1.d0
a=1.d0
v0=1.d0
rc=0.1d0
q=0.9d0
T=2.*!dpi
dt=T/10000
nstep=10*T/dt

;----------------------------------------------------;
; Vektorit lopputuloksille
;----------------------------------------------------;
r_vec=findgen(nstep)*0.d0
r_vec[0]=r0
dr_vec=r_vec
dr_vec[0]=dr0
z_vec=r_vec
z_vec[0]=z0
dz_vec=r_vec
dz_vec[0]=dz0
phi_vec=r_vec
phi_vec[0]=phi0
dphi_vec=r_vec
dphi_vec[0]=dphi0
timet=findgen(nstep)*dt

;------------------------------------------------------;
; RK4-looppi
;------------------------------------------------------;
time=0
rr=[r0,z0,phi0,dr0,dz0,dphi0]
for i=0,nstep-1 do begin
  der=deriv_cyl(time,rr)
  res=rk4(rr,der,time,dt,'deriv_cyl',/double)
  rr=res

;Tallennetaan vektoreihin
r_vec[i]=rr[0]
z_vec[i]=rr[1]
phi_vec[i]=rr[2]
dr_vec[i]=rr[3]
dz_vec[i]=rr[4]
dphi_vec[i]=rr[5]
endfor


pot=0.5d0*v0^2*alog(rc^2+r_vec^2+z_vec^2/q^2)
E=pot+0.5d0*(dr_vec^2+dz_vec^2+(r_vec*dphi_vec)^2)
L=r_vec^2*dphi_vec
x_vec=r_vec*cos(phi_vec)
y_vec=r_vec*sin(phi_vec)

!p.multi=[0,2,2]
nwin
plot,x_vec,y_vec,title=a
plot,timet,E,title='Energia  (RK4)',xr=[0,60]
plot,timet,L,title='Kulmaliikem‰‰r‰  (RK4)',xr=[0,60],yr=[0.99999999,1.00000001]

;Energian ja kulmaliikem‰‰r‰n suhteellinen muutos
remove,0,E
E_max=max(abs(E))
E_min=min(abs(E))

dE=(E_max-E_min)/E_max
print,'Energian suhteellinen muutos (RK4):',dE
print,'E_max=',E_max
print,'E_min=',E_min

remove,0,L
L_max=max(abs(L))
L_min=min(abs(L))

dL=(L_max-L_min)/L_max
print,'Kulmaliikem‰‰r‰n suhteellinen muutos (RK4):',dL
print,'L_max',L_max
print,'L_min',L_min
end
