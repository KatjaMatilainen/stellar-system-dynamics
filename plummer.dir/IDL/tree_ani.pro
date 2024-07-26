pro tree_ani,runid,fac0=fac0,step0=step0,nonew=nonew,full=full,cen=cen,$
linux=linux,nogas=nogas,maxstep=maxstep


if(n_params() le 0) then begin
    print,'tree_ani,runid'
    return
endif


fac=1.
if(keyword_set(fac0)) then fac=fac0


anifile=runid+'/fort.29'

close,1
swap_e=0
if(keyword_set(linux)) then swap_e=1


openr,1,/f77_unf,anifile,swap_e=swap_e

versio=0.
;readu,1,versio
;print,'versio=',versio

mui=0

step=0l
time=0.
npart=0l
ngas=0l
nass=0l
niraf0=0l
biraf=0.
ngas_use=ngas
if(keyword_set(nogas)) then ngas_use=0

on_ioerror,endi

while not eof(1) do begin

readu,1,step,time,npart,niraf0,biraf
readu,1,step,time,ngas,niraf0,biraf
readu,1,step,time,nass,niraf0,biraf

readu,1,xh,yh,zh,vxh,vyh,vzh

if(mui eq 0) then begin
    dens=lonarr(niraf0,niraf0)
    densz=dens
    dens_g=dens
    densz_g=densz
    dens_a=dens
    densz_a=densz
    niraf=1.*niraf0*fac
    if(not keyword_set(nonew)) then begin
        if(ngas_use le 0) then nwin,xs=2*niraf,ys=niraf+40
        if(ngas_use gt 0) then nwin,xs=2*niraf,ys=2*niraf+40
        endif
    xyouts,niraf,niraf+10,runid,/dev
    zero=lonarr(niraf/2,50)*0
endif
mui=mui+1

if(npart gt 20) then readu,1,dens
if(ngas  gt 20) then readu,1,dens_g
if(nass  gt 20) then readu,1,dens_a

if(npart gt 20) then readu,1,densz
if(ngas  gt 20) then readu,1,densz_g
if(nass  gt 20) then readu,1,densz_a

if(keyword_set(step0)) then begin
    if(mui mod step0 ne 0) then goto,skip
endif




if(ngas_use le 0) then begin
if(fac eq 1) then begin
    tvscl,alog(dens>.01),0,0
    tvscl,alog(densz>.01),niraf,0
endif

if(fac gt 1) then begin
    dens_u=congrid(dens,niraf,niraf)
    densz_u=congrid(densz,niraf,niraf)
    tvscl,alog(dens_u>.01),0,0
    tvscl,alog(densz_u>.01),niraf,0
endif

tv,zero,0,niraf
xyouts,10,niraf+10,'T= '+string(time,'(f6.1)'),/dev
if(keyword_set(cen)) then plots,niraf/2,niraf/2,/dev,psym=1,syms=2
endif


if(ngas_use gt 0) then begin
if(fac eq 1) then begin
    loadct,3,/sil
    tvscl,alog(dens>.01),0,niraf
    tvscl,alog(densz>.01),niraf,niraf
    loadct,9,/sil
    tvscl,alog(dens_g>.01),0,0
    tvscl,alog(densz_g>.01),niraf,0
endif

if(fac gt 1) then begin
    dens_u=congrid(dens,niraf,niraf)
    densz_u=congrid(densz,niraf,niraf)
    dens_gu=congrid(dens_g,niraf,niraf)
    densz_gu=congrid(densz_g,niraf,niraf)
    tvscl,alog(dens_u>.01),0,niraf
    tvscl,alog(densz_u>.01),niraf,niraf
    tvscl,alog(dens_gu>.01),0,0
    tvscl,alog(densz_gu>.01),niraf,0
endif

tv,zero,0,2*niraf
xyouts,10,2*niraf+10,'T= '+string(time,'(f6.1)'),/dev
if(keyword_set(cen)) then plots,niraf/2,niraf/2,/dev,psym=1,syms=2
endif



if(keyword_set(maxstep)) then begin
if(step ge maxstep) then goto,endi
endif

skip:

endwhile


endi:
;stop
end

