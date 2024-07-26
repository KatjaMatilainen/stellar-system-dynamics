;harj4_galaxies

;--------------------------------------------------------------
;JAKAUMA ra,dec

psopen,'harj4_galaxies.ps',/color

restore,'galaxies.save',/verbose 
nwin
plot,gara,gadec,xtitle='rektaskensio',ytitle='deklinaatio',psym=3

restore,'open_clusters.save',/verbose
oplot,ocra,ocdec,psym=6,col=2

plot_stamp,'harj4_galaxies'
psclose

;--------------------------------------------------------------
psopen,'harj4_galaxies_b.ps',/color
;JAKAUMA l,b

year=2000.
glactc,gara,gadec,year,gal,gab,1
plot,gal,gab,xtitle='galaktinen longitudi',ytitle='galaktinen latitudi',$
psym=3

glactc,ocra,ocdec,year,ocl,ocb,1
restore,'open_clusters.save',/verbose
oplot,ocl,ocb,psym=6,col=2

plot_stamp,'harj4_galaxies_b'


;--------------------------------------------------------------
psopen,'harj4_galaxies_c.ps',/color

mpc=vhelio/75.

ind=where(gadec gt 26.5 and gadec lt 32.5 and mpc gt 0)


plot,gara(ind),mpc(ind),xtitle='rektaskensio',ytitle='etaisyys [MPC]',xr=[10,18],xs=1,yr=[1,150],ys=1,title='26.5<dek<32.5',psym=6,syms=.2

plot_stamp,'harj4_galaxies_c.ps'
psclose


end
