;Ohjelma, joka plottaa rotaatiok‰yr‰n eksponentiaaliselle
;kiekolle, jonka M ja Re ovat tiedossa.
;K‰ytt‰‰ k‰yr‰n laskemiseen aliohjelmaa expdisc.pro


;R-vektori ja p‰‰tetyt vakiot:
R=findgen(500)*0.01d0
sig0=1.d0
Re=1.d0
Mtot=2.*!dpi*sig0*Re^2

;--------------------------------------------------------------------------;
;   Vertaillaan eksponentiaalista kiekkoa ja Toomre-Kuzminia
;--------------------------------------------------------------------------;

;Eksponentiaalinen kiekko:
expdisc,R,M=Mtot,Re=Re,vcirc=vcirc_exp

;Toomre-Kuzmin:
toomre,R,M=Mtot,Re=Re,vcirc=vcirc_toom

!p.multi=[0,2,2]
nwin
plot,R,vcirc_exp,title='Eksponentiaalinen kiekko',xtitle='R',ytitle='v_circ'
plot,R,vcirc_toom,title='Toomre-Kuzmin -kiekko',xtitle='R',ytitle='v_circ',col=2

;Molemmat samassa kuvassa:
plot,R,vcirc_exp,title='Vertailukuva',xtitle='R',ytitle='v_circ'
oplot,R,vcirc_toom,col=2

;---------------------------------------------------------------------------;
;Vertaillaan eksponentiaalista ja pallosymmetrist‰ jakaumaa
;---------------------------------------------------------------------------;
;Edit: siirretty omaksi tiedostokseen harj2_3.pro

end
