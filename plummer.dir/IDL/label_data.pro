pro label_data,x0,y0,labels,colors=colors,lines=lines,size=size, $
    psym=psym,symsize=symsize,t3d=t3d,title=title,t_color=t_color, $
    thick=thick,y_increment=y_inc,zvalue=zvalue,symbols=symbols,$
    center=center,box=box,len=len,shift_symbol=shift_symbol,hiplaus=hiplaus
;+
; NAME:		LABEL_DATA
;
; PURPOSE:	Label lines and plotting symbols on a plot.
;		The default is to draw lines of different linestyles
;		and colors, putting a label next to each line.
; CATEGORY:
; CALLING SEQUENCE:
;		label_data,xpos,ypos,labels
; INPUTS:
;		xpos=position of upper left-hand corner of label block
;		    along x axis, in normalized coordinates.
;		ypos=position of upper left-hand corner of label block
;		    along y axis, in normalized coordinates.
;		labels=n-element string array of labels.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;		colors=n-element array of color specifications for labels.
;		lines=n-element array of linestyle specifications for labels.
;		size=scalar specifiying the size of the type.
;		psym=n-element array of psym specifications for labels.
;		symsize=n-element array of symsize spec's for labels.
;		symbols=array of 'symbol' specifiers: each element of
;		    psym which is equal to 8 (user-defined symbol)
;		    must have a corresponding value for 'symbol' to be 
;		    used by the procedure SYMBOLS.  
;		    Examples:	psym=[8,8,8,8],symbols=[1,2,20,30]
;				psym=[1,2,8,8],symbols=[1,2]
;		title=scalar string for title of label block.
;		t_color=scalar specifying the color of the title.
;		thick=n-element array specifying the thick spec's for labels
;		y_increment=scalar specifying the space between consecutive
;		    labels, in percent of plot window.
;		center = set to center the title.
;		Graphics Keywords: t3d,zvalue
;		
; OUTPUTS:
; OPTIONAL OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;		D. L. Windt, AT&T Bell Laboratories, November 1989
;-
on_error,2

n_words=n_elements(labels)	; Get number of words
if n_params() lt 3 then begin
  print,'label_data: no words to label!'
  return
  endif

p_flag=0.
if n_elements(lines) ne n_words then $ ; If lines not specified...
  if n_elements(psym) ne n_words then $ ; and if psym not spec'd...
  lines=findgen(n_words) $      ; use default linestyles.
else begin 
    lines=fltarr(n_words)       ; else just plot symbols.
    p_flag=1.
endelse

;print,lines
if n_elements(colors) ne n_words then colors=fltarr(n_words)+1.
if n_elements(psym) ne n_words then psym=fltarr(n_words)
if n_elements(symsize) ne n_words then symsize=fltarr(n_words)+1.
if keyword_set(size) eq 0 then size=!p.charsize>1
if n_elements(thick) eq 0 then thick=!p.thick+fltarr(n_words)
if keyword_set(t3d) eq 0 then t3d=!p.t3d
if keyword_set(y_inc) eq 0 then y_inc=.075
if keyword_set(zvalue) eq 0 then zvalue=!z.crange(1)

x_len=!x.crange(1)-!x.crange(0)	; Get length of plot area
y_len=!y.crange(1)-!y.crange(0) ; Get height of plot area

x=fltarr(3)

l15=.15
if(keyword_set(len)) then l15=len

if p_flag then begin
    x(0)=(x0+l15)*x_len+!x.crange(0)	; Get position for plotting symbols.
    x(1)=x(0)
    endif else begin
	x(0)=x0*x_len+!x.crange(0)	; Get position for left end of lines.
	x(1)=x(0)+l15*x_len		; Get position for right end of lines.
	endelse

x(2)=x(1)+.025*x_len		; Get position for start of words.

if !x.type then x=10^x

y=y0*y_len+!y.crange(0)		; Get position for first word.

if keyword_set(title) then begin
    if keyword_set(t_color) eq 0 then t_color=!p.color 
    ytitle_pos=y+y_inc*y_len
    if !y.type then ytitle_pos=10^ytitle_pos
    if keyword_set(center) then align=.5 else align=0.
    xyouts,x(2),ytitle_pos,title, $
        size=size*1.05,t3d=t3d,z=zvalue,color=t_color,alignment=align
    endif

i=0
sindex=0
dypos=0
if(keyword_set(shift_symbol)) then dypos=shift_symbol

while i lt n_words do begin
  ypos=y-i*y_inc*y_len
  if !y.type then ypos=10^ypos
  
  ;plot symbols
  psym_use=0
  if psym(i) ne 0. then begin
      psym_use=psym(i)
      if(psym_use ge 8) then psym_use=8
      if(psym_use le -8) then psym_use=-8      
      if(not keyword_set(hiplaus)) then begin
          plots,[x(0),x(1)],[ypos,ypos],color=colors(i),thick=thick(i), $
            t3d=t3d,[zvalue,zvalue],psym=abs(psym_use),symsize=symsize(i)
      endif
      if(keyword_set(hiplaus) and lines(i) lt 0) then begin
          plots,[0.5*(x(0)+x(1))],[ypos],color=colors(i),thick=thick(i), $
            t3d=t3d,[zvalue,zvalue],psym=abs(psym_use),symsize=symsize(i)
      endif

  endif

  ;connect with lines
  if(lines(i) ge 0) then begin
      plots,[x(0),x(1)],[ypos,ypos],lines=lines(i),color=colors(i), $
        thick=thick(i),t3d=t3d,[zvalue,zvalue],psym=-abs(psym_use),symsize=symsize(i)
;      PRINT,'LINES',LINES(I),psym_use,dypos
 endif

  xyouts,x(2),ypos-dypos,labels(i),size=size,color=colors(i),t3d=t3d,z=zvalue
  i=i+1
  endwhile


 if(keyword_set(box)) then begin
 if(n_elements(box) eq 2) then begin
     x00=(x(0)-!x.crange(0))/x_len-.015
     y00=y0+.035
     x1=x00+box(0)
     y1=y00-box(1)
 endif
 if(n_elements(box) eq 4) then begin
     x00=box(0)
     y00=box(1)
     x1=box(2)
     y1=box(3)
 endif

 plots,[x00,x00,x1,x1,x00],[y00,y1,y1,y00,y00],lines=0,/normal
endif

return
end



