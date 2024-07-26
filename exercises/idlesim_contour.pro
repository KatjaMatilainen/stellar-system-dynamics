;*******************************************************************
;idlesim_contour.pro
;esimerkki väreillä täytetystä contour-kuvasta
;*******************************************************************

x=findgen(50)/50.*10 
data=cos(x)#(sin(x)*exp(-.25*x))  


clev=[-0.1, 0., 0.1, 0.3, 0.5, 0.7, 0.9]     ;contour-arvot\\
cin =[      2,  3,   4,   5,   6,   7,   8]  ;vastaavat värit\\
clab=[1,    1,  1,   0,   1,   0,   1]       ;label vai ei\\

;HUOM: ensimmaisen contour arvon alittavat alueet varilla 0 (musta ruudulla)

tek_color
contour,data,x,x,xrange=[0,10],yrange=[0,10],xstyle=1,ystyle=1,$ 
levels=clev,c_label=clab,c_col=cin,/fill

;piirretaaan paalle contour-kayrat valkoisella + arvot

contour,data,x,x,xrange=[0,10],yrange=[0,10],xstyle=1,ystyle=1,$
levels=clev,/over,c_label=clab,col=1,chars=2


end 
;*******************************************************************
