function deriv_h7,t,yy
common h7_common

n=n_elements(yy)/4.

r=yy[0l:n-1]
f=yy[n:2*n-1]
rp=yy[2*n:3*n-1]
fp=yy[3*n:4*n-1]

drdf=yy*0

apu_r=-r/(r^2+a^2)^1.5-2.d0*eps*a*r*(r^4-a^4)/(r^2+a^2)^4*cos(2*f)
rpp=apu_r+(omeb+fp)^2*r-2*lambda*rp

apu_f=-2.d0*eps*a*r/(r^2+a^2)^2*sin(2*f)
fpp=1.d0/r*(apu_f-2*rp*(omeb+fp))

drdf=[rp,fp,rpp,fpp]

return,drdf

end
