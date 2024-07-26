;Ohjelma, joka vertaa kesken‰‰n eksponentiaalisen kiekon ja plummer-spheren
;rotaatiok‰yri‰, kun niiden kumulatiivinen massa M(r) on sama

R=1.d-6*1.01d0^dindgen(2000)
;R=findgen(100)*0.01d0
sig0=1.d0
Re=1.d0
G=1.d0

;Typer‰sti nimetty ohjelma, mutta laskee exp-discille kumulatiivisen massan
;annetusta s‰teest‰ ja keskustiheydest‰.
vrot_exp,R,sig0=sig0,Re=Re,M=M

;T‰m‰ ohjelma laskee plummerille rotaatiok‰yr‰n annetusta massasta.
vrot_plummer,R,G=G,M=M,vcirc=vcirc_plum

;T‰m‰ ohjelma laskee exp-discille rotaatiok‰yr‰n
expdisc,R,M=M,Re=Re,vcirc=vcirc_exp

;Plotataan k‰yr‰t

!p.multi=[0,2,2]
nwin
plot,R,vcirc_exp,title='Eksponentiaalinen kiekko',xtitle='R',ytitle='v_circ'
plot,R,vcirc_plum,title='Plummer-pallo',xtitle='R',ytitle='v_circ',col=2

;Molemmat samassa kuvassa:
plot,R,vcirc_exp,xr=[0,6],yr=[0,2],title='Vertailukuva',xtitle='R',$
ytitle='v_circ'
oplot,R,vcirc_plum,col=2

end
