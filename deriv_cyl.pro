function deriv_cyl,t,rr
common h4_common,G,M,a,v0,rc,q
  r=rr[0]
  z=rr[1]
  phi=rr[2]
  dr=rr[3]
  dz=rr[4]
  dphi=rr[5]

apu=0.5d0*v0^2/(rc^2+r^2+z^2/q^2)

ddr=-apu*2*r+r*dphi^2
ddz=-apu*2*z/q^2
ddphi=-2./r*dphi*dr

der=rr*0

der[0]=dr
der[1]=dz
der[2]=dphi
der[3]=ddr
der[4]=ddz
der[5]=ddphi

return,der

end

