;*********************************************
pro box,x0,y0,x1,y1,color
polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],col=color
end
;*********************************************

;*********************************
pro colortable,close=close,tek=tek
;*********************************

;get window-index of active window 
old_window=!d.window
;print,old_window

; destroy colortable-window if /close
if(keyword_set(close)) then begin
wdelete,5
return
endif

;creat colortable-window
window,5,xsize=350,ysize=350

;not /TEK -> display all 256 colors

if(not keyword_set(tek)) then begin
zero=replicate(' ',20)
plot,[0,0],[0,0],/nodata,xtickname=zero,ytickname=zero,$
title='current colortable',xtit='',ytit='',$
xstyle=1,ystyle=1,xr=[5,345],yr=[15,345]
for i=0,15 do begin
for j=0,15 do begin
col=16*i+j
x1=10+i*20
x2=x1+10
y1=340-j*20
y2=y1-10

box,x1,y1,x2,y2,col
endfor
endfor

;/TEk -> display first 32 colors

endif else begin
zero=replicate(' ',20)
plot,[0,0],[0,0],/nodata,xtickname=zero,ytickname=zero,$
title='colors 0-31',xtit='',ytit='',$
xstyle=1,ystyle=1,xr=[5,345],yr=[15,345]
for i=0,3 do begin
for j=0,7 do begin
col=8*i+j
x1=20.+i*80
x2=x1+20
y1=330-j*40
y2=y1-20

box,x1,y1,x2,y2,col
xyouts,/data,x2+2,y2+4,string(col,'(i2)')
endfor
endfor
endelse


; make old window back to active window
; (if it did exist)

if(old_window ge 0) then wset,old_window
title=''

;return
end









