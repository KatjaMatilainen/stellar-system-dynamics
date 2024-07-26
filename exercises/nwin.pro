;*********************************************************
; nwin.pro
; opens new window 
; if device ne PS and if first plot in page
;*********************************************************

pro nwin,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,pix=pix,large=large,$
         small=small,huge=huge,info=info

if(keyword_set(info)) then begin
    print,'pro nwin,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos'
    print,'         ,pix=pix,large=large,small=small,huge=huge,info=info'
    print,'like IDL window, put works with PS-file and !p.multi'
    print,'HS ca 1990'
return
endif


if not keyword_set(xsize) then xsize=640
if not keyword_set(ysize) then ysize=512

if(keyword_set(large)) then begin
    xsize=640*1.6
    ysize=512*1.5
endif

if(keyword_set(huge)) then begin
    xsize=640*1.6*1.2
    ysize=512*1.5*1.2
endif

if(keyword_set(small)) then begin
    xsize=640*0.6
    ysize=512*0.6
endif


if(not keyword_set(pix)) then begin
    if(!d.name ne 'PS' and !d.name ne 'TEK' and !p.multi(0) eq 0) then  begin
        if (keyword_set(xpos) and keyword_set(ypos)) then $
          window,/free,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,retain=2
        if (not keyword_set(xpos) or not keyword_set(ypos)) then $
          window,/free,xsize=xsize,ysize=ysize,retain=2
    endif
endif


if(keyword_set(pix)) then begin
    if(!d.name ne 'PS' and !d.name ne 'TEK' and !p.multi(0) eq 0) then  begin
        if (keyword_set(xpos) and keyword_set(ypos)) then $
          window,/free,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,$
          /pix,retain=2        
        if (not keyword_set(xpos) or not keyword_set(ypos)) then $
          window,/free,xsize=xsize,ysize=ysize,/pix,retain=2
    endif
endif




end



