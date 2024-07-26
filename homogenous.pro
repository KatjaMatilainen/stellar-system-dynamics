;------------------------------------------------------------;
;Aliohjelma, joka laskee annetulle säteelle potentiaalin ja
;ratanopeuden käyttäen homogeenisen pallon tiheysjakaumaa
;------------------------------------------------------------;

pro homogenous,r,rho=rho,mass=mass,grav=grav,pot=pot,vrot=vrot

a=1.d0
G=1.d0
Mtot=4./3*!dpi*a^3

;Tiheysjakauma
index=n_elements(r)
rho=r*0.d0

for i=0,index-1 do begin
  if r(i) ge a then rho(i)=0.d0 else rho(i)=3*Mtot/(4*!dpi*a^3)
endfor

;Massa
v=r*0.d0
mass=r*0.d0

for i=1,index-1 do begin
  v(i)=4./3*!dpi*r(i)^3
  mass(i)=mass(i-1)+0.5d0*(rho(i)+rho(i-1))*(v(i)-v(i-1))
endfor
;print,'Kokonaismassa',max(mass)

;Gravitaatiovoima
grav=rho*0.d0
for i=1,index-1 do begin
  grav(i)=-G*mass(i)/r(i)^2
endfor

;Potentiaali
pot=rho*0.d0
for i=index-2,0,-1 do begin
  pot(i)=pot(i+1)+grav(i)*(r(i+1)-r(i))
endfor

;Nopeus ympyräradalla
vrot=rho*0.d0
for i=1,index-1 do begin
  vrot(i)=sqrt(r(i)*abs(grav(i)))
endfor

end
