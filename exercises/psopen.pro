;***************************************
;  PSOPEN.PRO
;  directs output to ps-file

   pro psopen,file,eps=eps,ctek=ctek,$
   color=color,vfont=vfont,sbox=sbox,copy=copy,up=up,right=right,$
   help=help,info=info,background=background

;**************************************************************

if(n_params() le 0 or keyword_set(info) or keyword_set(help)) then begin
print,' pro psopen,file,eps=eps,ctek=ctek,color=color,vfont=vfont'
print,' sbox=sbox,copy=copy,up=up,right=right,help=help,info=info'
print,' '
print,' direct output to file in landscape mode'
print,' /eps -> encapsulated ps :asks for more info'
print,' /color -> color ps'
print,' /ctek -> good at ctek-printer ?'
print,' /vfont -> use vector-fonts'
print,' /sbox -> force square shaped plot'
print,' /copy -> copy current colortable stretch (def=no)'
print,' up-> move upward by UP cm '
print,' right-> move right by RIGHT cm'
print,' /info or /help -> this message'
print,' HSalo 7.8.96  /appl1/heikki/UNAMIDL/psopen.pro'
print,' '
print,' background = color'
return
endif 
;**************************************************************

icopy=0
if(keyword_set(copy)) then icopy=1
set_plot,'ps',copy=icopy

device,FILENAME=file,/LANDSCAPE,bits_per_pixel=8,$
  color=0,scale=1.,encapsulated=0

; good values for sparc printer
xoff=1
yoff=24
if(keyword_set(up)) then xoff=xoff+up
if(keyword_set(right)) then yoff=yoff-right

device,xsize=22.5,ysize=18,yoff=yoff,xoff=xoff

if(keyword_set(ctek)) then device,yoff=27.
if(keyword_set(color)) then device,/color
sbox,0
if(keyword_set(sbox)) then sbox,1


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

!P.FONT=0
if(keyword_set(vfont)) then !p.font=-1

if(keyword_set(background)) then begin
    polyfill,[0,0,1,1,0],[0,1,1,0,0],col=background,/normal
endif

end







