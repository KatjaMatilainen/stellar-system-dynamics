function deriv_nsym,t,yy
common h5_common
  x=yy[0]
  y=yy[1]
  vx=yy[2]
  vy=yy[3]

  dydx=yy*0

ax=-v0^2*x/(rc^2+x^2+y^2/q^2)
ay=-v0^2*y/q^2/(rc^2+x^2+y^2/q^2)

dydx[0]=vx
dydx[1]=vy
dydx[2]=ax
dydx[3]=ay

return,dydx

end
