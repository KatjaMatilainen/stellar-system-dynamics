;*******************************************
; triad_read.pro
; reads triad-file containing asteroid data
; saves to triad.save
; HS 1992
;*******************************************

close,10
openr,10,'triad'


;***************************************************************
; triad contains following Asteroid data
; number a e sini wbar node dwbar dnode nres tomars tojup family

; number = asteroid number
; a      = semimajor axis
; e      = eccentricity
; sini   = sine of inclination
; wbar   = longitude of perihelion
; node   = longitude of ascending node
; dwbar  = precession rate of perihelion
; dnode  = precession rate of node
; nres   = index for secular resonance (if any)
; tomars = closest distance to mars
; tojup  = closest distance to jupiter
; family = family identification number (if any)

; elements are proper elements, corrected for periodic perturbations
; due to major planets
; distances are in astronomical units (AU)
; angles in degrees
; zero values mean no data: no family member
;			    no resonance
;                           proper eccentricity or inclination not determined

;****************************************************
; define arrays to hold elements

number=intarr(2000)
a=fltarr(2000)
e=a
sini=a
wbar=a
node=a
dwbar=a
dnode=a
nres=number
tomars=a
tojup=a
family=a

;*****************************************************
; read file line by line

lask=0
while not eof(10) do begin
line=''
readf,10,line
if(lask eq 35)then print,line
number(lask)=strmid(line,0,4)
a(lask)=float(strmid(line,9,5))

apu=strmid(line,15,4)              ;poimi e merkkijonona
e(lask)=float(strtrim(apu,2))   ;muuta liukuluvuksi -
                                ;strtrim(apu,2) poistaa ylimääräiset ' ':t
                                ;-> tyhjasta kentasta tulee ''
                                ;muutoin float antaa virheilmoituksen
     
apu=strmid(line,20,4)
sini(lask)=float(strtrim(apu,2))

apu=strmid(line,25,5)
wbar(lask)=float(strtrim(apu,2))

apu=strmid(line,31,5)
node(lask)=float(strtrim(apu,2))

apu=strmid(line,37,5)
dwbar(lask)=float(strtrim(apu,2))

apu=strmid(line,43,6)
dnode(lask)=float(strtrim(apu,2))

nres(lask)=long(strmid(line,50,2))

apu=strmid(line,53,6)
tomars(lask)=float(strtrim(apu,2))

apu=strmid(line,59,6)
tojup(lask)=float(strtrim(apu,2))

family(lask)=float(strmid(line,68,5))
lask=lask+1
endwhile

;***********************************************
; clean arrays by removing extra empty values

number=number(0:lask)
a=a(0:lask)
e=e(0:lask)
sini=sini(0:lask)
wbar=wbar(0:lask)
node=node(0:lask)
dwbar=dwbar(0:lask)
dnode=dnode(0:lask)
nres=nres(0:lask)
tomars=tomars(0:lask)
tojup=tojup(0:lask)
nres=nres(0:lask)
family=family(0:lask)

;*************************************************
; store to file triad.save
; can be read by restore,'triad.save'

save,file='triad.save',a,e,sini,wbar,node,$
dwbar,dnode,nres,tomars,tojup,family

print,'to restore triad-data'
print,'restore,''triad.save''


ind=where(a ne 0 and e ne 0 and sini ne 0 and wbar ne 0 and node ne 0, count)
close,1
openw,1,'triad_simple.dat'
printf,1,'ASTEROIDIEN RATAELEMENTTEJA (kts triad_read.pro)'
printf,1,'NUMBER   A    E   sinI    WBAR   NODE'
for i=0,count-1 do begin
    ii=ind(i)
    printf,1,ii,a(ii),e(ii),sini(ii),wbar(ii),node(ii),$
      f='(i6,f12.5,2f10.4,2f12.2)'
endfor
close,1
print,'simplified data written to "triad_simple.dat"'

;-------------------------------------------------------------------

nwin,xs=900,ys=400
histo_f,a,1,4,.01,/plot,psym=10

;resonances

period=[2.,3.,4.,5.]
arat=period^(-2./3.)
ajup=5.2
ares=ajup*arat

for i=0,n_elements(ares)-1 do begin
oplot,ares(i)*[1,1],[0,3],col=2
endfor


period=[3./2.,5./2., 7./2,9./2.]
arat=period^(-2./3.)
ajup=5.2
ares=ajup*arat

for i=0,n_elements(ares)-1 do begin
oplot,ares(i)*[1,1],[0,2.5],col=3,thick=2
endfor



period=[7./3.,8./3.]
arat=period^(-2./3.)
ajup=5.2
ares=ajup*arat
for i=0,n_elements(ares)-1 do begin
oplot,ares(i)*[1,1],[0,2],col=5,thick=2
endfor

end




