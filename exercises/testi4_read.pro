;***********************
; testi4_read.pro
; lukee tiedoston testi4
;*********************** 

print,'luetaan seuraavan nakoinen tiedosto'
spawn,'cat testi4'			

;spawn-komennolla voidaan antaa UNIX-komentoja IDL-proseduurista'
;vastaa interaktiivisen tilan $-merkkia (sitaei voi kayttaa ohjelmassa!)

print,' '
print,' '
print,'eli muotoa'
print,'x1,y1,z1	 string1'
print,'......'
print,'xn,yn,zn	 stringn'
print,' '


openr,1,'testi4'
x=fltarr(100) & y=x & z=x	;varauduttu sataan pisteeseen
laatu=STRARR(100)		;luodaan merkkitieto-taulukko
lask=0				
dums=''

while NOT EOF(1) DO BEGIN	;jatketaan jos ei ole kohdattu
				;tiedoston loppua
readf,1,dum1,dum2,dum3,dums	;luetaan skalaareina, silla vektorin
				;alkiota ei voi suoraan lukea
x(lask)=dum1
y(lask)=dum2
z(lask)=dum3
laatu(lask)=dums
lask=lask+1
endwhile

x=x(0:lask-1)
y=y(0:lask-1)
z=z(0:lask-1)
laatu=laatu(0:lask-1)
close,1

print,'TARKISTUS:'
print,'x'
print,x
print,'y'
print,y
print,'z'
print,z
print,'laatu'
print,laatu

end


