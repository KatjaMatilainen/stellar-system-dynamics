;***********************
; testi2_read.pro
; lukee tiedoston testi2
;*********************** 

print,'luetaan seuraavan nakoinen tiedosto'
spawn,'cat testi2'

print,' '
print,' '
print,'eli muotoa'
print,' x1  y1  z1'
print,' .....'
print,' xn  yn  zn'
print,' '

print,'luetaan tavalla 1: rivi kerrallaan'
print,'silla ei tiedeta rivien maaraa'
print,' '
openr,1,'testi2'
x=fltarr(100) & y=x & z=x	;varauduttu sataan pisteeseen
lask=0				
while NOT EOF(1) DO BEGIN	;jatketaan jos ei ole kohdattu
				;tiedoston loppua
readf,1,dum1,dum2,dum3		;luetaan skalaareina, silla vektorin
				;alkiota ei voi suoraan lukea
x(lask)=dum1
y(lask)=dum2
z(lask)=dum3
lask=lask+1
endwhile

x=x(0:lask-1)
y=y(0:lask-1)
z=z(0:lask-1)
close,1

print,'TARKISTUS:'
print,'x'
print,x
print,'y'
print,y
print,'z'
print,z

print,' '
print,'luetaan tavalla 2: vektoreina'
print,'oletetaan etta tiedetaan  tiedostossa olevan 4 rivia'
print,'kullakin 3 lukua'
print,' '
openr,1,'testi2'
n=4
temp=fltarr(3,n)	;luodaan valiaikainen taulukko
readf,1,temp		;luetaan siihen tiedostosta kaikki kerralla
x=temp(0,*)
y=temp(1,*)
z=temp(2,*)
close,1

print,'TARKISTUS:'
print,'x'
print,x
print,'y'
print,y
print,'z'
print,z

print,'huomaa etta ens tapa antoi vaakavektoreina jalkimmainen'
print,'pystyvektoreina (voidaan muuttaa TRANSPOSE-komennolla'

end








