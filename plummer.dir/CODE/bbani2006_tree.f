c**************************************************************************

	subroutine anitul(step,tnow,nparttot,pos)

c**************************************************************************
ccc	include 'sinc.for'



	integer step,nparttot,ngastot,nass
	real tnow
        common /hs/noutani,niraf,biraf 


	real pos(56864,3)
	real xx(56864),yy(56864),zz(56864)

c	generates max 600x600 tables for animation 
c	maximum, actual size determined by niraf
c	display area specified by biraf


c	DENSXY  	=	density projected to xy plane
c	DENSXZ  	=	density projected to xz plane

c	output into units 29

c***************************************************************
	
	integer dens(600,600)
	integer deng(600,600)
	integer dena(600,600)


c***************************************************************

c       define dummies

	nass=0
	ngastot=0
	time=tnow
	xh=0.
        yh=0.
	zh=0.
	vxh=0.
        vyh=0.
	vzh=0.

	do 1 i=1,nparttot
	   xx(i)=pos(i,1)
	   yy(i)=pos(i,2)
	   zz(i)=pos(i,3)
 1	continue


	if(niraf.lt.0) return
	if(biraf.lt.0) return
	if(niraf.eq.0) niraf=200
	if(biraf.eq.0) biraf=10.
        if(niraf.gt.600) niraf=600

	xlow=-biraf
	xhigh=biraf
	ylow=-biraf
	yhigh=biraf
	zlow=-biraf
	zhigh=biraf

	iorigin=1
	
	xorigo=xh
	yorigo=yh
	zorigo=zh
	vxorigo=vxh
	vyorigo=vyh
	vzorigo=vzh

	xwidth=xhigh-xlow
	ywidth=yhigh-ylow
	zwidth=zhigh-zlow

	write(29) step,time,nparttot,niraf,biraf
	write(29) step,time,ngastot,niraf,biraf
	write(29) step,time,nass,niraf,biraf

	i=1
	write(29) xh,yh,zh,vxh,vyh,vzh

c**********************************************************
c 	xy-projection: 
c**********************************************************

	do 10 i=1,niraf
	do 11 j=1,niraf
	   dens(i,j)=0
	   deng(i,j)=0
	   dena(i,j)=0
 11	continue
 10	continue

	niraf2=niraf/2
	if(nparttot.gt.20) then
	do 20 i=1,nparttot
	ix=(xx(i)-xorigo)/xwidth*niraf+niraf2
	iy=(yy(i)-yorigo)/ywidth*niraf+niraf2
	if(ix.gt.niraf.or.ix.lt.1) goto 20
	if(iy.gt.niraf.or.iy.lt.1) goto 20
	dens(ix,iy)=dens(ix,iy)+1
 20	continue	
	write(29) ((dens(i,j),i=1,niraf),j=1,niraf)
	endif

c	if(ngastot.gt.20) then
c	   do 21 i=1,ngastot
c	      ix=(xg(i)-xorigo)/xwidth*niraf+niraf2
c	      iy=(yg(i)-yorigo)/ywidth*niraf+niraf2
c	      if(ix.gt.niraf.or.ix.lt.1) goto 21
c	      if(iy.gt.niraf.or.iy.lt.1) goto 21
c	      deng(ix,iy)=deng(ix,iy)+1
c 21	   continue	
c	   write(29) ((deng(i,j),i=1,niraf),j=1,niraf)
c	endif
	
c	if(nass.gt.20) then
c	   do 22 i=1,nass
c	      ix=(xa(i)-xorigo)/xwidth*niraf+niraf2
c	      iy=(ya(i)-yorigo)/ywidth*niraf+niraf2
c	      if(ix.gt.niraf.or.ix.lt.1) goto 22
c	      if(iy.gt.niraf.or.iy.lt.1) goto 22
c	      dena(ix,iy)=dena(ix,iy)+1
c 22	   continue	
c	   write(29) ((dena(i,j),i=1,niraf),j=1,niraf)
c	endif

c       2d-no vertical projection

c        goto 999

c**********************************************************
c 	xz-projection: unit 29
c**********************************************************

	do 110 i=1,niraf
	   do 111 j=1,niraf
	      dens(i,j)=0
	      deng(i,j)=0
	      dena(i,j)=0
 111	   continue
 110	continue

	niraf2=niraf/2

	if(nparttot.gt.20) then
	   do 120 i=1,nparttot
	      ix=(xx(i)-xorigo)/xwidth*niraf+niraf2
	      iz=(zz(i)-zorigo)/zwidth*niraf+niraf2
	      if(ix.gt.niraf.or.ix.lt.1) goto 120
	      if(iz.gt.niraf.or.iz.lt.1) goto 120
	      dens(ix,iz)=dens(ix,iz)+1
 120	   continue	
	   write(29) ((dens(i,j),i=1,niraf),j=1,niraf)
	endif
	
c	if(ngastot.gt.20) then
c	   do 121 i=1,ngastot
c	      ix=(xg(i)-xorigo)/xwidth*niraf+niraf2
c	      iz=(zg(i)-zorigo)/zwidth*niraf+niraf2
c	      if(ix.gt.niraf.or.ix.lt.1) goto 121
c	      if(iz.gt.niraf.or.iz.lt.1) goto 121
c	      deng(ix,iz)=deng(ix,iz)+1
c 121	   continue	
c	   write(29) ((deng(i,j),i=1,niraf),j=1,niraf)
c	endif
	
c	if(nass.gt.20) then
c	   do 122 i=1,nass
c	      ix=(xa(i)-xorigo)/xwidth*niraf+niraf2
c	      iz=(za(i)-zorigo)/zwidth*niraf+niraf2
c	      if(ix.gt.niraf.or.ix.lt.1) goto 122
c	      if(iz.gt.niraf.or.iz.lt.1) goto 122
c	dena(ix,iz)=dena(ix,iz)+1
c 122	continue	
c	write(29) ((dena(i,j),i=1,niraf),j=1,niraf)
c	endif

 999	continue
	return
	end
