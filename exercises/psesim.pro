;***************************************
;  PSOPEN.PRO
;  directs output to ps-file
   pro psopen,file,color=color,vfont=vfont
;***************************************

if(n_params() eq 0) then begin
print,' pro psopen,file,color=color,vfont=vfont'
print,' directs output to ps-file FILE'
print,' landscape mode, 8-bits'
print,' /COLOR  -> create color-ps'
print,' /VFONT  -> use vector-fonts (default=hardware'
print,' use DEVICE-command to override settings'
return
endif

set_plot,'ps'

device,FILENAME=file,/LANDSCAPE,bits_per_pixel=8,color=0,scale=1.
device,xsize=22.5,ysize=18,yoff=24,xoff=1	;sparcilla hyva

if(keyword_set(color)) then device,/color

!P.FONT=0
if(keyword_set(vfont)) then !p.font=-1
return
end



;***************************************
;  PSCLOSE.PRO
;  closes ps-file
   pro psclose,junk=junk
;***************************************

device,/CLOSE_FILE
!P.FONT=-1

if(keyword_set(junk)) then begin
spawn,'lpr junk'
print,'sent to printer'
endif
set_plot,'X'
end


;********************************************
; sbox.pro
; creates square-sized plot on default window
; with aspect ratio 1.25
;********************************************

pro sbox,size

if n_params() eq 0 then begin
print,' pro sbox,size
print,' square-plot on default window
print,' size determines plot size '
print,' size=0 -> restores default settings'
return
endif

xwidth=0.64*size*1.08939
ywidth=0.8*size*1.08939
x1=.5-xwidth/2.
x2=.5+xwidth/2.
y1=.5-ywidth/2.
y2=.5+ywidth/2.

!p.position=[x1,y1,x2,y2]
if(size eq 0) then !p.position=0

return
end

