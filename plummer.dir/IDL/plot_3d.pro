pro plot_3d,x,y,z,wid0=wid0,name=name,syms0=syms0

if(n_params() le 0) then begin
   print,'plot_3d,x,y,z,wid0=wid0,name=name,syms0=syms0'
   return
endif

wid=max([abs(x),abs(y),abs(z)])
if(keyword_set(wid0)) then wid=wid0
tek_color
TVLCT,R,G,B,/GET


syms=0.01*wid
if(keyword_set(syms0)) then syms=syms0*wid


oOrb = OBJ_NEW('orb', COLOR=[0, 0, 255])  
oOrb->Scale, syms, syms, syms  
oSymbol = OBJ_NEW('IDLgrSymbol', oOrb)  


XPLOT3D, X, Y, Z, NAME=name, $  
   SYMBOL=oSymbol, THICK=2 ,lines=6

end
