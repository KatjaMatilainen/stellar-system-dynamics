pro use_xplot3d,x,y,z,wid0=wid0

if(n_params() le 0) then begin
print,'use_xplot3d,x,y,z,wid=wid'
print,' '
print,'use xplot3d to make an interactive 3D plot of'
print,'orbits defined by x(*,i),y(*,i),z(*,i) i=0,npart-1'
print,'wid -> define limits of plotting cube from -lim to lim'
return
endif

npart=n_elements(x(0,*))

wid=max([abs(x),abs(y),abs(z)])
if(keyword_set(wid0)) then wid=wid0
tek_color
TVLCT,R,G,B,/GET

for i=0,npart-1 do begin
col=[R(I+2),G(I+2),B(I+2)]
if(i eq 0) then xplot3d,x(*,i),y(*,i),z(*,i),col=col,$
                xr=[-1,1]*wid,yr=[-1,1]*wid,zr=[-1,1]*wid

if(i ne 0) then xplot3d,x(*,i),y(*,i),z(*,i),col=col,/over
endfor
return
end
