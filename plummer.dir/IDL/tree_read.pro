;********************************************************
; tree_read.pro  <-- TREAD.PRO 4.12.90 - 28.3.95
; reads projection files (formatted) created by tree.exe 
; stores particle positions and halo-positions (dummies)
; into save- and saveindex-files
; to be processed with same routines as bbcode-output
;********************************************************

pro tree_read,runid,runid1,min,max,dire=dire,$
              restore=restore,noplot=noplot,ind1=ind1,ind2=ind2,auto=auto

  if(n_params() le 0) then begin
     print,'pro tread,runid,runid1,min,max,restore=restore'
     print,'read treecode data in files runidNNN (e.g. TEST/SNAP)'
     print,'(from SNAPmin - SNAPmax)'
     print,'output pos/vel  to idl-savefiles runid1sxyzMMM  (MMM=NNN-1)'
     print,'       mds      to runid1mmxyz'
     print,'mean values etc to runid1sxyzindex'
     print,' '
     print,'/restore-> restore from runid1sxyzNNN-files'
     print,'ind1,ind2 -> store and plot separately'
     print,'/auto -> all runid/SNAPNNN -> runid/runidsxyzMMM'
     print,'noplot=0 -> display vx,vy,vz disp'
     print,'            90% 75% 50% 25 % 10% of mass + 100%/10'
     print,'            initial and final r,vx diagram'
     print,' '
     print,"Example:  tree_read,'plummer_example2010',/auto"
     print,'TJD2010  HS'
     return
  endif




  if(keyword_set(auto)) then begin
     min=1
     max=1000
     runid1=runid+'/'+runid
     runid=runid+'/'+'SNAP'
  endif

  if(keyword_set(dire)) then runid=dire+'/'+runid
  if(keyword_set(dire)) then runid1=dire+'/'+runid1


  ngal=2
  mhalob=fltarr(ngal)
  mdisk=mhalob
  mgas=mhalob
  xhalo=fltarr(100,ngal)
  yhalo=xhalo
  zhalo=xhalo
  vxhalo=xhalo
  vyhalo=xhalo
  vzhalo=zhalo

  xh=fltarr(ngal)
  yh=xh
  zh=xh
  vxh=xh
  vyh=xh
  vzh=xh

  nna=lonarr(ngal)
  nnb=nna

  savtable=strarr(100)
  timetable=fltarr(100)
  colltable=lonarr(100)
  steptable=lonarr(100)
  vxxtable=fltarr(100)
  vyytable=fltarr(100)
  vzztable=fltarr(100)
  mrrtable=fltarr(100)
  drrtable=fltarr(100)
  cxtable=fltarr(100)
  cytable=fltarr(100)
  cztable=fltarr(100)
  f10table=fltarr(100)
  f25table=fltarr(100)
  f50table=fltarr(100)
  f75table=fltarr(100)
  f90table=fltarr(100)
  f100table=fltarr(100)

  savindex=runid1+'sxyzindex'
  massfile=runid1+'mmxyz'

  for ifile=min-1,max-1 do begin

     iname=ifile+1
     if(iname le 999) then file=runid+''+string(iname,'(i3)')
     if(iname le 99) then file=runid+'0'+string(iname,'(i2)')
     if(iname le 9) then file=runid+'00'+string(iname,'(i1)')
     
     if(ifile le 999) then sfile=runid1+'sxyz'+string(ifile,'(i3)')
     if(ifile le 99) then sfile=runid1+'sxyz'+string(ifile,'(i2)')
     if(ifile le 9) then sfile=runid1+'sxyz'+string(ifile,'(i1)')
     

;************************************************
;  RESTORING FILES
;************************************************

     if(keyword_set(restore)) then begin
        
        close,1
        openr,1,file,err=stat
        if(stat ne 0) then goto,skip_reading2
        close,1
        
        print,'reading from ... ',sfile
        restore,sfile
        print,time
        goto,skip_reading
     endif

;************************************************
;  READING FILES
;************************************************


     close,1
     openr,1,file,err=stat
     if(stat ne 0) then goto,skip_reading2
     
     print,'reading from ... ',file
     readf,1,ntul
     readf,1,ndim
     readf,1,time
     print,time,ntul
     
     temp=fltarr(1,ntul)
     readf,1,temp               ;massat
     if(ifile eq min-1) then begin
        mds=temp
        save,filename=massfile,mds
     endif

     
     temp=fltarr(3,ntul)
     readf,1,temp               ;paikat
     xx=temp(0,*)
     yy=temp(1,*)
     zz=temp(2,*)
     
     readf,1,temp               ;nopeudet
     vx=temp(0,*)
     vy=temp(1,*)
     vz=temp(2,*)
     
     nna(0)=0
     nna(1)=ntul-1
     nnb(0)=ntul-1
     nnb(1)=ntul-1
     
;tarpeetonta lukea
;temp=fltarr(1,ntul)
;readf,1,temp	;eps
;hlp,temp
     
     xx=reform(xx)
     yy=reform(yy)
     zz=reform(zz)
     vx=reform(vx)
     vy=reform(vy)
     vz=reform(vz)
     ncoll=0
     
     step=0
     cx=mean(xx)
     cy=mean(yy)
     cz=mean(zz)
     dist=sqrt((xx-cx)^2+(yy-cy)^2+(zz-cz)^2)
     
     print,sfile
     save,time,step,ncoll,xh,yh,zh,vxh,vyh,vzh,xx,yy,zz,vx,vy,vz,ngal,nna,nnb,$
          FILENAME=sfile,dist,cx,cy,cz
     
skip_reading:
     
     
     timetable(ifile)=time
     savtable(ifile)=sfile
     
     print,' '
     print,'-------------------------'
     print,'           POSITIONS:                            VELOCITIES'
     PRINT,'       MEAN         STDEV                    MEAN        STDEV'
     print,'X ',mean(xx),stdev(xx),'   ',mean(vx),stdev(vx)
     print,'Y ',mean(yy),stdev(yy),'   ',mean(vy),stdev(vy)
     print,'Z ',mean(zz),stdev(zz),'   ',mean(vz),stdev(vz)
     print,'-------------------------'
     
     
     cx=mean(xx)
     cy=mean(yy)
     cz=mean(zz)
     dist=sqrt((xx-cx)^2+(yy-cy)^2+(zz-cz)^2)
     order=sort(dist)
     ntul=n_elements(xx)
     f10=dist(order(ntul*0.10))
     f25=dist(order(ntul*0.25))
     f50=dist(order(ntul*0.50))
     f75=dist(order(ntul*0.75))
     f90=dist(order(ntul*0.90))
     f100=max(dist)
     
     vxxtable(ifile)=stdev(vx)
     vyytable(ifile)=stdev(vy)
     vzztable(ifile)=stdev(vz)
     cxtable(ifile)=mean(xx)
     cytable(ifile)=mean(yy)
     cztable(ifile)=mean(zz)
     mrrtable(ifile)=mean(dist)
     drrtable(ifile)=stdev(dist)

     f10table(ifile)=f10
     f25table(ifile)=f25
     f50table(ifile)=f50
     f75table(ifile)=f75
     f90table(ifile)=f90
     f100table(ifile)=f100
     
     if(ifile eq 0) then begin
        vx0=vx
        dist0=dist
     endif
     
  endfor

skip_reading2:
  max=ifile-1

  nparttot=ntul
  ngastot=0
  nindex=ifile-1


;if(max(timetable) eq 0) then timetable=findgen(100)*0.2
  timetable=timetable(0:nindex)
  savtable=savtable(0:nindex)
  vxxtable=vxxtable(0:nindex)
  vyytable=vyytable(0:nindex)
  vzztable=vzztable(0:nindex)
  cxtable=cxtable(0:nindex)
  cytable=cytable(0:nindex)
  cztable=cztable(0:nindex)
  mrrtable=mrrtable(0:nindex)
  drrtable=drrtable(0:nindex)
  f10table=f10table(0:nindex)
  f25table=f25table(0:nindex)
  f50table=f50table(0:nindex)
  f75table=f75table(0:nindex)
  f90table=f90table(0:nindex)
  f100table=f100table(0:nindex)

  save,runid,savtable,nindex,ngal,mhalob,mdisk,mgas,nparttot,ngastot,$
       xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo,steptable,timetable,colltable,$
       FILENAME=savindex,vxxtable,vyytable,vzztable,cxtable,cytable,cztable,$
       mrrtable,drrtable

;****
;plot

  if(not keyword_set(noplot)) then begin
     !p.multi=[0,2,2]
     nwin
     plot,timetable,vxxtable,psym=-4,col=1,title=runid1,xtit='time',ytit='vdisp'
     oplot,timetable,vyytable,psym=-2,col=2,lines=2
     oplot,timetable,vzztable,psym=-1,col=3,lines=1
     
;plot,timetable,mrrtable,psym=-4,col=1,title=runid1,xtit='time',ytit='rdist'
;oplot,timetable,drrtable,psym=-2,col=2,lines=2
     plot,timetable,f90table,psym=0,col=2,lines=1,$
          title=runid1,xtit='time',ytit='R(mf)'
     oplot,timetable,f100table/10,psym=-4,col=1
     oplot,timetable,f75table,psym=0,lines=1,col=3
     oplot,timetable,f50table,psym=0,lines=1,col=4
     oplot,timetable,f25table,psym=0,lines=1,col=5
     oplot,timetable,f10table,psym=0,lines=1,col=6
     

     dmax=max([abs(dist0),abs(dist)])
     vmax=max([abs(vx0),abs(vx)])
     !x.range=[0,max(dist)]
     !y.range=[-vmax,vmax]
     plot,dist0,vx0,psym=3,col=1,title=runid1+' INIT',xtit='dist',ytit='vx'
     plot,dist,vx,psym=3,col=2,title=runid1+'FINAL',xtit='dist',ytit='vx'
  endif

;****************************************
; IND-files
;****************************************

  for round=1,2 do begin
     
     if(round eq 1) then begin
        if(not keyword_set(ind1)) then goto,skip_ind
     endif
     
     if(round eq 2) then begin
        if(not keyword_set(ind2)) then goto,skip_ind
     endif
     
     savtable=strarr(100)
     timetable=fltarr(100)
     colltable=lonarr(100)
     steptable=lonarr(100)
     vxxtable=fltarr(100)
     vyytable=fltarr(100)
     vzztable=fltarr(100)
     mrrtable=fltarr(100)
     drrtable=fltarr(100)
     cxtable=fltarr(100)
     cytable=fltarr(100)
     cztable=fltarr(100)
     f10table=fltarr(100)
     f25table=fltarr(100)
     f50table=fltarr(100)
     f75table=fltarr(100)
     f90table=fltarr(100)
     f100table=fltarr(100)

     for ifile=min-1,max-1 do begin
        if(ifile le 999) then sfile=runid1+'sxyz'+string(ifile,'(i3)')
        if(ifile le 99) then sfile=runid1+'sxyz'+string(ifile,'(i2)')
        if(ifile le 9) then sfile=runid1+'sxyz'+string(ifile,'(i1)')
        restore,sfile
        
        if(round eq 1) then begin
           xx=xx(ind1)
           yy=yy(ind1)
           zz=zz(ind1)
           vx=vx(ind1)
           vy=vy(ind1)
           vz=vz(ind1)
        endif
        if(round eq 2) then begin
           xx=xx(ind2)
           yy=yy(ind2)
           zz=zz(ind2)
           vx=vx(ind2)
           vy=vy(ind2)
           vz=vz(ind2)
        endif

        timetable(ifile)=time
        savtable(ifile)=sfile
        cx=mean(xx)
        cy=mean(yy)
        cz=mean(zz)
        dist=sqrt((xx-cx)^2+(yy-cy)^2+(zz-cz)^2)
        order=sort(dist)
        ntul=n_elements(xx)
        f10=dist(order(ntul*0.10))
        f25=dist(order(ntul*0.25))
        f50=dist(order(ntul*0.50))
        f75=dist(order(ntul*0.75))
        f90=dist(order(ntul*0.90))
        f100=max(dist)
        
        vxxtable(ifile)=stdev(vx)
        vyytable(ifile)=stdev(vy)
        vzztable(ifile)=stdev(vz)
        cxtable(ifile)=mean(xx)
        cytable(ifile)=mean(yy)
        cztable(ifile)=mean(zz)
        mrrtable(ifile)=mean(dist)
        drrtable(ifile)=stdev(dist)
        
        f10table(ifile)=f10
        f25table(ifile)=f25
        f50table(ifile)=f50
        f75table(ifile)=f75
        f90table(ifile)=f90
        f100table(ifile)=f100

        if(ifile eq 0) then begin
           vx0=vx
           dist0=dist
        endif
     endfor
     nparttot=ntul
     ngastot=0
     nindex=ifile-1

     timetable=timetable(0:nindex)
     savtable=savtable(0:nindex)
     vxxtable=vxxtable(0:nindex)
     vyytable=vyytable(0:nindex)
     vzztable=vzztable(0:nindex)
     cxtable=cxtable(0:nindex)
     cytable=cytable(0:nindex)
     cztable=cztable(0:nindex)
     mrrtable=mrrtable(0:nindex)
     drrtable=drrtable(0:nindex)
     f10table=f10table(0:nindex)
     f25table=f25table(0:nindex)
     f50table=f50table(0:nindex)
     f75table=f75table(0:nindex)
     f90table=f90table(0:nindex)
     f100table=f100table(0:nindex)
     
     if(round eq 1) then begin
        title=runid1+'IND: '+string(min(ind1),'(I5)')+string(max(ind1),'(I5)')
        save,runid,savtable,nindex,ngal,mhalob,mdisk,mgas,nparttot,ngastot,$
             xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo,$
             steptable,timetable,colltable,$
             FILENAME=savindex+'1',vxxtable,vyytable,$
             vzztable,cxtable,cytable,cztable,$
             mrrtable,drrtable,ind1,title
     endif
     
     if(round eq 2) then begin
        title=runid1+'IND: '+string(min(ind2),'(I5)')+string(max(ind2),'(I5)')
        save,runid,savtable,nindex,ngal,mhalob,mdisk,mgas,nparttot,ngastot,$
             xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo,$
             steptable,timetable,colltable,$
             FILENAME=savindex+'2',vxxtable,vyytable,vzztable,$
             cxtable,cytable,cztable,$
             mrrtable,drrtable,ind2,title
     endif
     
     
     if(not keyword_set(noplot)) then begin
        defplot
        !p.multi=[0,2,2]
        nwin
        plot,timetable,vxxtable,psym=-4,col=1,$
             title=title,xtit='time',ytit='vdisp'
        oplot,timetable,vyytable,psym=-2,col=2,lines=2
        oplot,timetable,vzztable,psym=-1,col=3,lines=1
        
        plot,timetable,f90table,psym=0,col=2,$
             lines=1,title=runid1,xtit='time',ytit='R(mf)'
        oplot,timetable,f100table/10,psym=-4,col=1
        oplot,timetable,f75table,psym=0,lines=1,col=3
        oplot,timetable,f50table,psym=0,lines=1,col=4
        oplot,timetable,f25table,psym=0,lines=1,col=5
        oplot,timetable,f10table,psym=0,lines=1,col=6
        
        dmax=max([abs(dist0),abs(dist)])
        vmax=max([abs(vx0),abs(vx)])
        !x.range=[0,max(dist)]
        !y.range=[-vmax,vmax]
        plot,dist0,vx0,psym=3,col=1,title=title+' INIT',xtit='dist',ytit='vx'
        plot,dist,vx,psym=3,col=2,title=title+'FINAL',xtit='dist',ytit='vx'
     endif
     
skip_ind:
  endfor                        ;round

  defplot

end
