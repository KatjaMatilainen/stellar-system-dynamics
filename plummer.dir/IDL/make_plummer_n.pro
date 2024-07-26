;*********************************************************
;  make_plummer_n
;  Create radius and velocity vectors from the
;  distribution function corresponding to a Plummer sphere
;  HS 28.03.07  TJD2007
;*********************************************************


      G=1.  
for iround=1,2 do begin   


   if(iround eq 1) then begin
      npart=40l
      ngal=npart
      b=100.
      M=40.
   endif

   if(iround eq 2) then begin
      npart=500l
      b=1.
      M=1.
      epsi=0.05*b
   endif

;unifromly distributed random numbers
  s1=randomu(seed,npart)
  s2=randomu(seed,npart)
  s3=randomu(seed,npart)

;******************************************************************
;particle positions: 
;   absolute values according to Plummer distribution
;   directions isotropic
;******************************************************************
  timer,/start
  rad=b*sqrt(s1^(2./3.)/(1.-s1^(2./3.)))

  fii=!pi*2.*s2
  sinf=sin(fii)
  cosf=cos(fii)

  cost=s3*2.-1.
  sint=sqrt(1.-cost^2)
  
  xx=rad*sint*cosf
  yy=rad*sint*sinf
  zz=rad*cost

  timer,/stop,dt
  print,'POSITIONS required ',dt

;--------------------------------------------------
;check the density profile by comparing to analytic
;--------------------------------------------------
;tabulation --> rhotest(rtest)
;               lkmtest(rtest) = number 
;               -> 1/sqrt(number) gives an error estimate

  rtest=findgen(40)/20.*b
  rhotest=rtest*0.
  drtest=rtest(1)-rtest(0)
  lkmtest=rtest*0.+1.

; theoretical Plummer distribution

  rhotheory= 3.*m/4./!pi/b^3*(1.+(rtest/b)^2)^(-5./2.)

  for i=0,n_elements(rtest)-2 do begin
      ind=where(rad ge rtest(i) and rad lt rtest(i+1),count)
      vol=4.*!pi/3.*(rtest(i+1)^3-rtest(i)^3)
      if(count ge 1) then rhotest(i)=count/vol*m/npart    
      if(count ge 1) then lkmtest(i)=count        
  endfor
  rel_err=1./sqrt(lkmtest)

  nwin
  plot,(rtest+drtest*.5)/b,rhotest,$
    xtitle='r/b',ytitle='rho',yr=[0,1.5*rhotheory(0)],psym=10
  oplot,(rtest+drtest*.5)/b,rhotest,psym=10,col=3
  
  for i=0,n_elements(rtest)-1 do begin
      oplot,(rtest(i)+drtest*.5)/b*[1,1],$
        rhotest(i)*[1.-rel_err(i),1.+rel_err(i)],col=3
  endfor
  oplot,rtest/b,rhotheory,col=2,lines=2
  label_data,0.5,0.9,['MC','Theory'],lines=[0,2],col=[3,2]


;****************************************************************
;velocities
;****************************************************************

  timer,/start
;construct G(theta)

  theta=(findgen(1000)+0.5)/1000.*!pi/2.
  sint8=sin(theta)^8
  sint10=sin(theta)^10
  ftheta=sint8-sint10
  gtheta=theta*0.
  apu=total(ftheta)
  summ=0.
  for i=999,0,-1 do begin
      summ=summ+ftheta(i)
      gtheta(i)=summ
  endfor
  gtheta=gtheta/apu
  
  s4=randomu(seed,npart)
  s5=randomu(seed,npart)
  s6=randomu(seed,npart)
  
;for each s4, search the index of the nearest tabulated value
  indt=lonarr(npart)
  for i=0l,npart-1 do begin
      dev=min(abs(s4(i)-gtheta))
      ind=where(abs(s4(i)-gtheta) eq dev)
      indt(i)=ind(0)
  endfor

  psi=G*M/b/sqrt(1.+rad^2/b^2)
  vrad=sqrt(2.*psi)*cos(theta(indt))

  fii=!pi*2.*s5
  sinf=sin(fii)
  cosf=cos(fii)
  cost=s6*2.-1.
  sint=sqrt(1.-cost^2)

  vx=vrad*sint*cosf
  vy=vrad*sint*sinf
  vz=vrad*cost

  enekin=0.5*(vx^2+vy^2+vz^2)
  ee=psi-enekin
  eptable=ee
;--------

  vv=sqrt(vx^2+vy^2+vz^2)
  rr=sqrt(xx^2+yy^2+zz^2)
  print,'max r,v',max(rr),max(vv)

  timer,/stop,dt
  print,'VELOCITIES required ',dt


print,'-------------------------'
print,'POSITIONS:'
PRINT,'  MEAN      STDEV'
print,'X ',mean(xx),stdev(xx)
print,'Y ',mean(yy),stdev(yy)
print,'Z ',mean(zz),stdev(zz)

print,'-------------------------'
print,'VELOCITIES:'
PRINT,'  MEAN      STDEV'
print,'X ',mean(vx),stdev(vx)
print,'Y ',mean(vy),stdev(vy)
print,'Z ',mean(vz),stdev(vz)
print,' '


if(iround eq 1) then begin
   xgal=xx
   ygal=yy
   zgal=zz
   vxgal=vx
   vygal=vy
   vzgal=vz
endif

if(iround eq 2) then begin
   mtotal=m
   pmass=xx*0.+mtotal/npart
   eps=xx*0.+epsi
endif
endfor



xx0=-999
yy0=-999.
zz0=-999.
vx0=-999
vy0=-999.
vz0=-999.
pmass0=-999.
eps0=-999.

for igal=0,ngal-1 do begin

   xx0=[xx0,xgal(igal)+xx]
   yy0=[yy0,ygal(igal)+yy]
   zz0=[zz0,zgal(igal)+zz]
   vx0=[vx0,vxgal(igal)+vx]
   vy0=[vy0,vygal(igal)+vy]
   vz0=[vz0,vzgal(igal)+vz]
   pmass0=[pmass0,pmass]
   eps0=[eps0,eps]
endfor

xx=xx0(1:*)
yy=yy0(1:*)
zz=zz0(1:*)
vx=vx0(1:*)
vy=vy0(1:*)
vz=vz0(1:*)
pmass=pmass0(1:*)
eps=eps0(1:*)
mtotal=total(pmass)

xc=total(pmass*xx)/mtotal
yc=total(pmass*yy)/mtotal
zc=total(pmass*zz)/mtotal

vxc=total(pmass*vx)/mtotal
vyc=total(pmass*vy)/mtotal
vzc=total(pmass*vz)/mtotal

xx=xx-xc
yy=yy-yc
zz=zz-zc
vx=vx-vxc
vy=vy-vyc
vz=vz-vzc
ndim=3
time=0.

npart=n_elements(xx)

TREE_FILE='plummer_example_n.TREEBI'

close,1
openw,1,tree_file
printf,1,npart
printf,1,ndim
printf,1,time
for i=0l,npart-1 do begin
    printf,1,pmass(i)
endfor
for i=0l,npart-1 do begin
    printf,1,xx(i),yy(i),zz(i)
endfor
for i=0l,npart-1 do begin
    printf,1,vx(i),vy(i),vz(i)
endfor
for i=0l,npart-1 do begin
    printf,1,eps(i)
endfor
close,1

endi:
end













