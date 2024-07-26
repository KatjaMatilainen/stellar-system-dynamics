;harj3_cassini

restore,'saturn.save'
loadct,0
nwin,xs=3*512,ys=512
tvscl,rebin(red,512,512),0
tvscl,rebin(green,512,512),1
tvscl,rebin(blue,512,512),2
xyouts,256+512*0,400,'RED',/dev,ali=0.5
xyouts,256+512*1,400,'GREEN',/dev,ali=0.5
xyouts,256+512*2,400,'BLUE',/dev,ali=0.5
plot_stamp,'harj5_cassini'




end
