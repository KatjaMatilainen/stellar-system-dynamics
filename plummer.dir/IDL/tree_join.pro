;**************************
;tree_join.pro
;**************************

pro tree_join,filein,fileout,xc0=xc0,yc0=yc0,zc0=zc0,$
              vxc0=vxc0,vyc0=vyc0,vzc0=vzc0,$
              xx=xx,yy=yy,zz=zz,vvx=vvx,vvy=vvy,vvz=vvz,pmass=pmass

  if(n_params() le 0) then begin
     print,'-----------------------------------------'
     print,' pro tree_join,filein,fileout'
     print,'-----------------------------------------'
     print,' join together TREEBI.files given in filein'
     print,' output to fileout (none-> no)'
     print,' '
     print,' relative shift:'
     print,' xc=xc,yc=yc,zc=zc,vxc=vxc,vyc=vyc,vzc=vzc'
     print,' created coordinates:  '
     print,' xx=xx,yy=yy,zz=zz,vvx=vvx,vvy=vvy,vvz=vvz,pmass=pmass'
     print,'-----------------------------------------'
     print,' EXAMPLEs:'
     print,"  infile = replicate('plummer_example2010.TREEBI',2)"
     print,"  outfile= 'circular_example.TREEBI'"
     print,"  xcir=5.  & vycir=sqrt(2./xcir)"
     print,' tree_join,infile,outfile,xc0=[0,xcir],vyc=[0,vycir]'
     print,'-----------------------------------------'
     print,"  infile = replicate('plummer_example2014.TREEBI',2)"
     print,"  outfile= 'plummer_example2014_collision.TREEBI'"
     print,' tree_join,infile,outfile,xc0=[0,10],vxc=[0,-2]'
     print,'-----------------------------------------'
     print,'HS 100510 /250314'
     return
  endif

  nfile=n_elements(filein)
  xc=fltarr(nfile)
  yc=xc
  zc=xc
  vxc=xc
  vyc=xc
  vzc=xc

  if(keyword_set(xc0)) then xc=xc0
  if(keyword_set(yc0)) then yc=yc0
  if(keyword_set(zc0)) then zc=zc0

  if(keyword_set(vxc0)) then vxc=vxc0
  if(keyword_set(vyc0)) then vyc=vyc0
  if(keyword_set(vzc0)) then vzc=vzc0


  for i=0,nfile-1 do begin

     file=filein(i)
     close,1
     openr,1,file

     print,'reading file...',file

     readf,1,npartfin
     readf,1,ndim
     readf,1,time
     xt=findgen(npartfin)
     yt=xt
     zt=xt
     vxt=xt
     vyt=xt
     vzt=xt
     pmasst=xt
     epst=xt

     for j=0l,npartfin-1 do begin
        readf,1,d1
        pmasst(j)=d1
     endfor
     for j=0l,npartfin-1 do begin
        readf,1,d1,d2,d3
        xt(j)=d1+xc(i)
        yt(j)=d2+yc(i)
        zt(j)=d3+zc(i)
     endfor

     for j=0l,npartfin-1 do begin
        readf,1,d1,d2,d3
        vxt(j)=d1+vxc(i)
        vyt(j)=d2+vyc(i)
        vzt(j)=d3+vzc(i)
     endfor

     for j=0l,npartfin-1 do begin
        readf,1,d1
        epst(j)=d1
     endfor
     close,1

     print,'-------------------------'
     print,'ifile=',i
     print,'file=',file
     print,'npart',npartfin
     print,'mass',total(pmasst)
     print,'xx',mean(xt),stdev(xt)
     print,'yy',mean(yt),stdev(yt)
     print,'zz',mean(zt),stdev(zt)
     print,'vx',mean(vxt),stdev(vxt)
     print,'vy',mean(vyt),stdev(vyt)
     print,'vz',mean(vzt),stdev(vzt)
     print,'-------------------------'

     if(i eq 0) then begin
        xx=xt
        yy=yt
        zz=zt
        vx=vxt
        vy=vyt
        vz=vzt
        pmass=pmasst
        eps=epst
     endif

     if(i gt 0) then begin
        xx=[xx,xt]
        yy=[yy,yt]
        zz=[zz,zt]
        vx=[vx,vxt]
        vy=[vy,vyt]
        vz=[vz,vzt]
        pmass=[pmass,pmasst]
        eps=[eps,epst]
     endif

  endfor                        ;i=ifile

;*************************
;mass-center
;*************************

  npart=n_elements(xx)

  cm=total(pmass)
  cx=total(pmass*xx)/cm
  cy=total(pmass*yy)/cm
  cz=total(pmass*zz)/cm
  cvx=total(pmass*vx)/cm
  cvy=total(pmass*vy)/cm
  cvz=total(pmass*vz)/cm

  xx=xx-cx
  yy=yy-cy
  zz=zz-cz
  vx=vx-cvx
  vy=vy-cvy
  vz=vz-cvz


  print,'-------------------------'
  print,'COMBINED'
  print,'npart',npart
  print,'mass',total(pmass)
  print,'xx',mean(xx),stdev(xx)
  print,'yy',mean(yy),stdev(yy)
  print,'zz',mean(zz),stdev(zz)
  print,'vx',mean(vx),stdev(vx)
  print,'vy',mean(vy),stdev(vy)
  print,'vz',mean(vz),stdev(vz)
  print,'-------------------------'



;*************************
;write to file
;*************************

  if(fileout ne 'none') then begin
     npartfin=npart
     tree_file=fileout

     print,'write to ',tree_file

     close,1
     openw,1,tree_file
     printf,1,long(npartfin)
     printf,1,long(ndim)
     printf,1,time
     for i=0l,npartfin-1 do begin
        printf,1,pmass(i)
     endfor
     for i=0l,npartfin-1 do begin
        printf,1,xx(i),yy(i),zz(i)
     endfor

     for i=0l,npartfin-1 do begin
        printf,1,vx(i),vy(i),vz(i)
     endfor

     for i=0l,npartfin-1 do begin
        printf,1,eps(i)
     endfor
     close,1
  endif

end
