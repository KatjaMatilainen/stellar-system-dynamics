;-------------------------------------------------------
pro siniplot,x,y,a,b,c,$
oplot0=oplot0,color0=color0,linestyle0=linestyle0,$
title0=title0
;-------------------------------------------------------

if(n_params() le 0) then begin
    print,'siniplot,x,y,a,b,c'
    print,'plottaa sinikayran y=A*sin(b*x+c)'
    print,'x maaritelty kutsuvassa paaohjelmassa, samoin a,b,c'
    print,'y palautetaan kutsuvaan ohjelmaan' 
    print,'keywordit:'
    print,'title        ->tee tasta plotin otsikko (default= a,b,c arvot'
    print,'color=color                 -> kayta varia color (def=1)'
    print,'linestyle=linestyle         -> viivatyyli (def=0)'
    print,'oplot=1  (sama kuin /oplot) -> plottaa entisen paalle'
    return
endif

y=a*sin(b*x+c)

col=1
if(keyword_set(color0)) then col=color0

line=0
if(keyword_set(linestyle0)) then line=linestyle0

ff='(f8.3)'
title='a,b,c='+string(a,ff)+string(b,ff)+string(c,ff)
if(keyword_set(title0)) then title=title0

if(not keyword_set(oplot0)) then begin
    plot,x,y,line=line,col=col,xtitle='x',ytitle='y=a*sin(b*x+c)',$
      title=title
endif

if(    keyword_set(oplot0)) then begin
    oplot,x,y,line=line,col=col
endif

end

;--------------------------------------------------------------
;PAAOHJELMA ALKAA TASTA
;--------------------------------------------------------------
;harj1.pro


;psopen,'harj1_siniplot.ps'   ;tulostus ps-tiedostoon

!p.multi=[0,1,2]          ;piirretaan 4 kuvaa

nwin                      ;avaa uuden ikkunan
                          ;HUOM: ei ole IDL:n oma komento
                           

x=findgen(1001)/1000.*2*!pi*10    
!x.range=[0,max(x)]
!x.style=1.


a=1 & b=1 & c=0
title='amplitudi muuttuu'
siniplot,x,y,a,b,c,title=title
siniplot,x,y,a*.50,b,c,/opl,col=2,lines=1
siniplot,x,y,a*.25,b,c,/opl,col=3,lines=2

a=1 & b=1 & c=0
title='jakso muuttuu'
siniplot,x,y,a,b,c,title=title
siniplot,x,y,a,b*.5,c,/opl,col=2,lines=1
siniplot,x,y,a,b*.25,c,/opl,col=3,lines=2

!p.multi=0
!x.style=0
!x.range=0


;plot_stamp,'harj1_siniplot.ps' ;kirjoitetaan ikkunaan infoa
;psclose         ;sulkee ps-tulostuksen
                ;HUOM: eivat ole IDL:n omia komentoja

end






