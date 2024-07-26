;**********************************************
   pro histo_f,x,x1,x2,dx,xx,yy,gg,$
   plot=plot,oplot=oplot,auto=auto,noscale=noscale,color=color,psym=psym,$
   gauss=gauss,surf=surf,dens=dens,cum=cum
;**********************************************
; distribution of x
; sampling area x1 -> x2 with bin dx
; return yy=f(xx) distribution normalized to unity
;**********************************************

if(n_params() eq 0) then begin
    print,'pro histo_f,x,x1,x2,dx,xx,yy,gg'
    print,'x=input values'
    print,'histogram from x1 to x2 with step dx'
    print,'xx,yy return calculated values'
    print,'gg returns corresponding gaussian fit (if /gauss)'
    print,'/plot      -> plot data'
    print,'/oplot     -> oplot data'
    print,'/auto      -> automatic scaling with 50 part/bin'
    print,'auto=nbins -> automatic scaling with nbin bins'
    print,'/noscale   -> do not scale distribution (default: normalized area)'
    print,'color=col  -> plot with given color'
    print,'psym =sym  -> plot with given symbol type'
    print,'surf=var   -> returns surface density
    print,'dens=var   -> returns volume density
    print,'cum=var    -> return  cumulative mass'
    print,'HSalo ca 1990'
return
end


nx=n_elements(x)

if(keyword_set(auto)) then begin
    x1=min(x)
    x2=max(x)
    if(auto gt 1) then  dx=(x2-x1)/auto
    if(auto le 1) then  dx=50*(x2-x1)/nx
endif

nbins=(x2-x1)/dx
xx=findgen(nbins)*dx+0.5*dx+x1	; centers of bins
num=lonarr(nbins)
num=histogram(x,min=x1,max=x2,binsize=dx)
yy=num(0:nbins-1)


if(not keyword_set(noscale)) then yy=1.*num/nx/dx

if(not keyword_set(color)) then color=!p.color
if(not keyword_set(psym)) then psym=!p.psym
if(keyword_set(plot) or keyword_set(auto)) then plot,xx,yy,col=color,psym=psym
if(keyword_set(oplot)) then oplot,xx,yy,col=color,psym=psym

if(keyword_set(gauss)) then begin
    sd=stdev(x)
    me=mean(x)
    gg=1./sqrt(2.*!pi)/sd*exp(-0.5*(xx-me)^2/sd^2)
    if(keyword_set(noscale)) then gg=gg*nx*dx
    oplot,xx,gg,col=3,psym=0,linestyle=2
    print,'mean =',me
    print,'stdev=',sd
endif


area2=(xx+0.5*(xx(1)-xx(0)))^2
area1=(xx-0.5*(xx(1)-xx(0)))^2
area=!pi*(area2-area1)
surf=yy/area

area2=(xx+0.5*(xx(1)-xx(0)))^3
area1=(xx-0.5*(xx(1)-xx(0)))^3
area=4.*!pi/3.*(area2-area1)
dens=yy/area

cum=xx*0
cums=0.
cum(0)=0.
for i=1l,n_elements(xx)-1 do begin
    cums=cums+yy(i)
    cum(i)=cums
endfor

return
end









