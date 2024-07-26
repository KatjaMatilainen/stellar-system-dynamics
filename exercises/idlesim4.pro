;----------------------------------------------------------------------
;idlesim4.pro
;esimerkki pinnan maarittelemisesta sisakkaisten for-silmukoiden avulla:
;sijoittamalla erilaisia funktioita saadaan mielivaltaisia z=z(x,y) pintoja
;  Huom! aiemmissa esimerkeissa on luotu 2-ulotteisia taulukoita 
;ulkotulo-operaation  \# avulla,
;jolloin z(x,y) on ollut aina muotoa z(x,y)=f(x)*f(y).
;----------------------------------------------------------------------

x=fltarr(51,51) & y=x   & z= x

FOR i=0,50 DO BEGIN 
    FOR j=0,50 DO BEGIN 
        x(i,j)=float(i) 
        y(i,j)=float(j) 
    ENDFOR
ENDFOR

x=x/5.-5.                       ;x on vektori, sisältää arvoja väliltä -5,5
y=y/5.-5.                       ;y samoin
xlab=findgen(51)/5.-5.          ;akselien merkintää varten
ylab=xlab 

z=exp(-0.125*(x^2+y^2))*sin(2*y)*cos(y)

loadct,3                        ;otetaan käyttöön eräs IDL-paleteista
shade_surf,z,xlab,ylab          ;plotataan varjostettuna pintana

end
;---------------------------------------------------------------------------
