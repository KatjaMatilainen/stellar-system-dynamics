npart=500
ndim=3
b=2.d0
M=1.d0
G=1.d0
eps=findgen(npart)*0.+0.01*b
time=0.

pmass=findgen(npart)*0.+M/npart

s1=randomu(seed,npart)
s2=randomu(seed,npart)
s3=randomu(seed,npart)
s4=randomu(seed,npart)
s5=randomu(seed,npart)
s6=randomu(seed,npart)

;Paikat pallokoord:
r=b*sqrt(s1^(2./3)/(1-s1^(2./3)))
f=2.*!dpi*s2
cth=2.*s3-1
th=acos(cth)

;Karteesisessa koord:
xx=r*sin(th)*cos(f)
yy=r*sin(th)*sin(f)
zz=r*cos(th)

;print,'x',xx
;print,'y',yy
;print,'z',zz

;Nopeudet:

vx=xx*0.
vy=yy*0.
vz=zz*0.

;Jotta saat nopeudet, laske ensin potentiaali.

pot=G*M/b*(1+r^2/b^2)^(-0.5)

;Nopeus:

;Tarvitaan apumuuttuja
apu=findgen(10^5)/(10^5)

v=sqrt(2.*pot)*cos(apu)

;Plottaa x,y -tasossa
nwin
plot,xx,yy,psym=3,xrange=[-10,10],yrange=[-10,10],iso=1

close,1
openw,1,'plummer_oma.TREEBI' 

 printf,1,npart                              ;N
 printf,1,ndim                               ;dimension = 3
 printf,1,time 
 for i=0l,npart-1 do begin
    printf,1,pmass(i)                       ;masses
 endfor
 for i=0l,npart-1 do begin
    printf,1,xx(i),yy(i),zz(i)              ;positions
 endfor
 for i=0l,npart-1 do begin
    printf,1,vx(i),vy(i),vz(i)              ;velocities
 endfor
 for i=0l,npart-1 do begin
    printf,1,eps(i)                         ;softening lengths
 endfor

close,1

end
