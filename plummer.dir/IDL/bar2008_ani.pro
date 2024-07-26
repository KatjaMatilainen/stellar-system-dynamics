pro bar2008_ani,runid,fac0=fac0,step0=step0,$
                nonew=nonew,full=full,cen=cen,omeb=omeb,linux=linux,z=z,$
                mpeg=mpeg,nframes0=nframes0,wait=wait,orbmax=orbmax,dire=dire

if(n_params() le 0) then begin
    print,'---------------------------------------------------------'
    print,' bar2008_ani,runid,fac0=fac0,step0=step0,$'
    print,'             nonew=nonew,full=full,cen=cen,omeb=omeb,linux=linux'
    print,'              circle=circle'
    print,'---------------------------------------------------------'
    print,' runid = simulation to show (expanded to runid/runid_ani)'
    print,' KEYWORDS:'
    print,'  fac    =  enlargement factor (def=1)'
    print,'  step   =  show every step stores (def=1)'
    print,'  /nonew -> no new window opened'
    print,'  /full  -> runid is not expanded'
    print,'  /cen   -> mark center'
    print,'  omeb   =  rotate with this angular speed (def=0)'
    print,'            -999 -> use omeb stored in animation-file (if it exists)'
    print,'  /linux -> read files made in other OS'
    print,'z=z,wait=wait'
    print,' '
    print,' HS 2006 / mpeg added 20.05.2008'
    print,'---------------------------------------------------------'
    return
endif

fac=1.
if(keyword_set(fac0)) then fac=fac0

anifile=runid+'/'+runid+'_ani'
if(keyword_set(full)) then anifile=runid+'_ani'
if(keyword_set(dire)) then anifile=dire+'/'+runid+'/'+runid+'_ani'


close,1

swap_e=0
if(keyword_set(linux)) then swap_e=1
openr,1,/f77_unf,anifile,swap_e=swap_e

versio=0.
readu,1,versio
print,'versio=',versio
if(versio eq 2) then begin
    readu,1,omebar_use
    print,omebar_use

    pbar=1.

    if(keyword_set(omeb)) then begin
        if(omeb le -999) then omeb=omebar_use
        pbar=2.*!pi/omeb
    endif
endif

mui=0

step=0l
time=0.
npart=0l
ngas=0l
nass=0l
niraf0=0l
biraf=0.

on_ioerror,endi

iframe=0
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
            nwin,xs=2*niraf,ys=niraf+40
        endif

        add=''
        if(keyword_set(z)) then add=' z-type'
        
        xyouts,niraf,niraf+10,runid+add,/dev
        zero=lonarr(niraf/2,50)*0
        print,niraf0,biraf


;mpeg-output
        if(keyword_set(mpeg)) then begin
            mpeg_file=runid+'.mpg'
            nframes=100
            if(keyword_set(nframes0)) then nframes=nframes0
            mpegid=mpeg_open([2*niraf,niraf+40])
        endif


    endif
    mui=mui+1

    if(npart gt 20) then readu,1,dens
    if(ngas  gt 20) then readu,1,dens_g
    if(nass  gt 20) then readu,1,dens_a
    
;bar output are 2D
    goto,eit
    if(npart gt 20) then readu,1,densz
    if(ngas  gt 20) then readu,1,densz_g
    if(nass  gt 20) then readu,1,densz_a
eit:
    
    if(keyword_set(step0)) then begin
        if(mui mod step0 ne 0) then goto,skip
    endif

    if(keyword_set(omeb)) then begin
        dens=rot(dens,omeb*time*!radeg,1.,niraf0/2-1.5,niraf0/2-1.5,/int)
        dens_g=rot(dens_g,omeb*time*!radeg,1.,niraf0/2-1.5,niraf0/2-1.5,/int)

        if(keyword_set(z)) then begin
            dens=reverse(dens)
            dens_g=reverse(dens_g)
        endif
    endif



    if(fac eq 1) then begin
        loadct,3,/sil
        tvscl,alog(dens>.01),0,0
        loadct,8,/sil
        tvscl,alog(dens_g>.01),niraf,0
    endif

    if(fac gt 1) then begin
        dens_u=congrid(dens,niraf,niraf)
        dens_gu=congrid(dens_g,niraf,niraf)
        loadct,3,/sil
        tvscl,alog(dens_u>.01),0,0
        loadct,8,/sil
        tvscl,alog(dens_gu>.01),niraf,0
    endif

    tv,zero,0,niraf
    xyouts,10,niraf+10,'T= '+string(time/pbar,'(f6.2)'),/dev
    if(keyword_set(cen)) then plots,niraf/2,niraf/2,/dev,psym=1,syms=2

 if(keyword_set(wait)) then wait,wait
;mpeg-output
        if(keyword_set(mpeg)) then begin
            windex=!d.window
            mpeg_put,mpegid,frame=iframe,/order,/color,window=windex
            if(iframe mod 20 eq 0) then print,'frame=',iframe

            if(iframe gt nframes) then goto,endi
        endif


if(keyword_set(orbmax)) then begin
if(time/pbar ge orbmax) then goto,endi
endif
iframe=iframe+1    


skip:
endwhile



endi:
;mpeg-output
        if(keyword_set(mpeg)) then begin
            mpeg_save,mpegid,file=mpeg_file
            mpeg_close,mpegid
            cmd='ls -altr '+mpeg_file
            spawn,cmd
        endif

end

