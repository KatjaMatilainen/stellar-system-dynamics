;***********************
; testi3_read.pro
; lukee tiedoston testi3
;*********************** 

print,'luetaan seuraavan nakoinen tiedosto'
spawn,'cat testi3'		

print,' '
print,' '
print,'eli muotoa'
print,'n'
print,'x1 x2 ... xn'
print,'y1 y2 ... yn'
print,'z1 z2 ... zn'
print,' '

openr,1,'testi3'
readf,1,n		;luetaan alkioiden lukumaara
x=fltarr(n)		;luodaan taulukkot
y=x
z=x

readf,1,x		;luetaan tiedostosta vektori x
readf,1,y		;luetaan tiedostosta vektori y
readf,1,z		;luetaan tiedostosta vektori z
close,1

print,' '
print,'TARKISTUS:'
print,'x'
print,x
print,'y'
print,y
print,'z'
print,z

end



