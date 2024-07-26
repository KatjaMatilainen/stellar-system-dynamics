;********************************
; PSDUMP.PRO
; dumps the active screen to ps-file
; 
;*******************************
pro psdump,file,savef=savef,white=white,portrait=portrait,tek=tek,eps=eps,$
inv=inv,wbg=wbg,color=color,nocheck=nocheck,hp=hp,scalef=scalef,$
bits=bits,a3=a3,ctek=ctek,xoff=xoff,yoff=yoff,wl=wl,kepler=kepler,tycho=tycho,$
copy=copy,help=help,info=info,hewlett=hewlett,up=up,right=right,temp=temp,$
full=full,true=true

if(keyword_set(help) or keyword_set(info)) then begin
print,'************************************************************'
print,'psdump      -> dump active screen to printer'
print,'psdump,file -> dump active screen psfile'
print,'
print,'default printer is sparc except if:' 
print,'/color -> dump to cjeta4'
print,'/a3    -> dump to cjeta3'
print,'/hp    -> turbops'
print,'/hewwlet -> hp in 3rd floor'
print,'/ctek  -> ctek'
print,' '
print,'/tek -> change tek-palette to bw'
print,'/inv -> invert colorscale'
print,'/white -> force white background (!d_ncoll -> white)'
print,' WL=value -> all values > WL into white'
print,' '
print,'savef=file -> dump directly and save to file'
print,'bits=value -> use bits_per_pixel=value (def=8 except for 1 for /tek)'
print,' '
print,'default scaling: certain to fit into page:'
print,'xsize<1.33*ysize -> scale=1'
print,'xsize>1.33*ysize -> scale=ysize/xsize* 1.33'
print,'scalef=value of scaling (if automatic scaling makes too small plots'
print,' '
print,'/portrait (def=landscape)'
print,'/eps -> encapsulated postricpt (asks for additional paramaters)'
print,'/copy  -> retain current colortable (def=no) NEW FEATURE !'

print,' '
print,'/appl1/heikki/UNAMIDL/psdump.pro ; HS 1990/14.6.96' 
print,'***************************************************************'
return
endif

;hp converted to turbops
;and lpq removed 1.9.92

old_index=!d.window
xsize0=!d.x_size
ysize0=!d.y_size
scale=1
if(n_params() eq 0) then begin
print,'dump directly to printer'
file='junks.ps'
if(keyword_set(color)) then file='junks.clps'
endif

if(1.*xsize0/ysize0 gt 1.33) then begin
scale=1.33/xsize0*ysize0
endif

if(keyword_set(a3)) then scale=scale*1.6
if(keyword_set(scalef)) then scale=scalef
print,scale

if(not keyword_set(true)) then temp=tvrd(0,0,xsize0,ysize0)
if(    keyword_set(true)) then begin
temp=tvrd(0,0,xsize0,ysize0,/true)
temp=color_quan(temp,1,rr,gg,bb,col=256)
tvlct,rr,gg,bb
endif

if(keyword_set(full)) then temp=bytscl(temp)

if(keyword_set(inv)) then temp=255b-temp

if(keyword_set(wbg)) then begin
max=max(temp)
temp(where(temp eq max))=1
temp(where(temp eq 0))=max
endif

if(keyword_set(white)) then begin
max=max(temp)
index=where(temp eq max)
temp(index)=255b
endif

if(keyword_set(wl)) then begin
index=where(temp ge wl)
temp(index)=255b
endif

if(keyword_set(tek)) then begin
indexb=where(temp ne 0b,countb) ;draw-region -> black
indexw=where(temp eq 0b,countw) ;background  -> white
if(countb ge 1) then temp(indexb)=0b
if(countw ge 1) then temp(indexw)=255b
endif

if(keyword_set(nocheck)) then begin
endif else begin
window,1,xsize=xsize0,ysize=ysize0,title=file
tv,temp
wshow,!d.window,0,/iconic
endelse

icopy=0
if(keyword_set(copy)) then icopy=1
psopen,file,copy=icopy
device,color=0
if(keyword_set(tek)) then device,bits_per_pixel=1
if(keyword_set(bits)) then device,bits_per_pixel=bits
if(keyword_set(color)) then device,/color
if(keyword_set(portrait)) then device,/PORTRAIT
if(keyword_set(color)) then begin
device,yoff=26.6

endif

if(keyword_set(xoff)) then device,xoff=0.75*2.54+xoff
if(keyword_set(yoff)) then device,yoff=27.+yoff
device,yoff=26.6

if(keyword_set(up)) then device,xoff=1+up
if(keyword_set(right)) then device,yoff=24-right

if(keyword_set(eps)) then begin
device,/ENCAPSULATED
print,'ENCAPSULATED POSTRIPT CHOSEN'
if(n_elements(eps) lt 5) then begin
print,'PORTRAIT-MODE (1-yes)'& read,vast & if(vast eq 1) then device,/PORTRAIT
print,'XSIZE (cm) ' & read,xsize & device,XSIZE=xsize
print,'YSIZE (cm) ' & read,ysize & device,YSIZE=ysize
print,'XOFF (cm) ' & read,xoff & device,XOFFSET=xoff
print,'YOFF (cm) ' & read,yoff & device,YOFFSET=yoff
endif
if(n_elements(eps) eq 5) then begin
if(eps(0) eq 1) then begin
device,/portrait
print,'PORTRAIT-MODE'
endif
if(eps(0) ne 1) then print,'LANDSCAPE-MODE'
device,xsize=eps(1)
device,ysize=eps(2)
device,xoffset=eps(3)
device,yoffset=eps(4)
print,'xsize   =',eps(1)
print,'ysize   =',eps(2)
print,'xoffset =',eps(3)
print,'xoffset =',eps(4)
endif

endif else begin
device,SCALE=scale
endelse

if(keyword_set(true)) then tvlct,rr,gg,bb
if(keyword_set(true)) then hlp,temp
tv,temp
psclose
wset,old_index

if(keyword_set(hewlett)) then begin
command='lpr -Php '+psfile+' &'
spawn,command
goto,endi
endif

if(keyword_set(savef)) then begin
print,'save image to file'
save,filename=file+'.psdump',temp
psfile=file+'.psdump'

if(keyword_set(kepler)) then begin
command='lpr -Plk '+psfile+' &'
spawn,command
endif

if(keyword_set(tycho)) then begin
command='lpr -Plt '+psfile+' &'
spawn,command
endif
endif



if not keyword_set(hp) then begin
if(file eq 'junks.ps') then begin
if(not keyword_set(kepler)) then command='lpr -h junks.ps &'
if(    keyword_set(kepler)) then command='lpr junks.ps &'
spawn,command
endif
endif else begin
if(file eq 'junks.ps') then begin
command='lpr -Pturbops junks.ps &'
spawn,command
;command='lpq -Pturbops'
;spawn,command
endif
endelse

if(file eq 'junks.clps') then begin
if(not keyword_set(a3)) then command='lpr -Pcjeta4 junks.clps &'
if(    keyword_set(a3)) then command='lpr -Pcjeta3 junks.clps &'
if(    keyword_set(ctek)) then command='lpr -Ptek junks.clps &'
spawn,command
if(    keyword_set(a3)) then spawn,'lpq -Pcjeta3'
if(not keyword_set(a3)) then spawn,'lpq -Pcjeta4'
if(    keyword_set(ctek)) then spawn,'lpq -Ptek'
endif
endi:
end














