;harj3_polyfit

nwin

x=(findgen(101)/100.-.5)*10
y=x^2+randomn(seed,101)
apu=poly_fit(x,y,2,yfit)
title='yfit='+string(apu(0))+'  +  '+string(apu(1))+' *x  +'+string(apu(2))+' *x^2'


plot,x,y,psym=6,syms=.5,title=title,xtitle='x',ytitle='y'
oplot,x,yfit,lines=2


print,title
plot_stamp,'harj3_polyfit'

end
