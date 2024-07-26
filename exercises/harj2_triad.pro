;-----------------------------------------
;harj2_triad.pro
;-----------------------------------------

;----------------------------------------------
;luetaan asteroidi-data
;----------------------------------------------
;kts myos triad_read.pro joka on lukenut alkuperaisen tiedoston

number=lonarr(2000)
a=fltarr(2000)
e=a
sini=a
wbar=a
node=a

close,10
openr,10,'triad_simple.dat'


;luetaan ensin kommenttirivit

line=''
readf,10,line & print,line
readf,10,line & print,line


lask=0
while not eof(10) do begin
    i=lask
    readf,10,number1,a1,e1,sini1,wbar1,node1
    number(i)=number1
    a(i)=a1
    e(i)=e1
    sini(i)=sini1
    wbar(i)=wbar1
    node(i)=node1
    lask=lask+1
endwhile

;siivotaan tyhjat pois
number=number(0:lask)
a=a(0:lask)
e=e(0:lask)
sini=sini(0:lask)
wbar=wbar(0:lask)
node=node(0:lask)



;----------------------------------------------
;leikitaan datalla
;----------------------------------------------

psopen,'harj2_family.ps'

nwin
!p.multi=[0,1,2]
plot,a,e,xtitle='keskietaisyys',ytitle='eksentrisyys',psym=3,xr=[2,4],$
title='Asteroidiperheet'

ind=where((a gt 3 and a lt 3.03) and (e gt 0.05 and e lt 0.1))
oplot,a(ind),e(ind),psym=3,col=2


plot,a,sini,xtitle='keskietaisyys',ytitle='sin(inklinaatio)',psym=3,xr=[2,4]
oplot,a(ind),sini(ind),psym=3,col=2
!p.multi=0

plot_stamp,'harj2_family.ps'
psclose















psopen,'harj2_kirkwood.ps'
nwin,xs=900,ys=400
!x.title='KESKIETAISYYS'
!y.title='LUKUMAARA'
!p.title='Kirkwoodin aukot'
histo_f,a,1,4,.01,/plot,psym=10
defplot
;resonances

period=[2.,3.,4.,5.]
s_period=['2:1','3:1','4:1','5:1']
arat=period^(-2./3.)
ajup=5.2028
ares=ajup*arat

for i=0,n_elements(ares)-1 do begin
    oplot,ares(i)*[1,1],[0,3],col=2
    xyouts,ares(i),3.1,ali=.5,s_period(i)
endfor


period=[3./2.,5./2., 7./2,9./2.]
s_period=['3:2','5:2','7:2','9:2']
arat=period^(-2./3.)
ajup=5.2028
ares=ajup*arat

for i=0,n_elements(ares)-1 do begin
oplot,ares(i)*[1,1],[0,2.5],col=3,thick=2
    xyouts,ares(i),2.6,ali=.5,s_period(i)
endfor



period=[7./3.,8./3.,10./3.]
s_period=['7:3','8:3','10:3']
arat=period^(-2./3.)
ajup=5.2028
ares=ajup*arat
for i=0,n_elements(ares)-1 do begin
oplot,ares(i)*[1,1],[0,2],col=5,thick=2
    xyouts,ares(i),2.1,ali=.5,s_period(i)
endfor

plot_stamp,'harj2_kirkwood.ps'
psclose
end
