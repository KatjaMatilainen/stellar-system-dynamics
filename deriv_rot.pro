function deriv_rot,t,yy
common h6_common

x=yy[0]
y=yy[1]
vx=yy[2]
vy=yy[3]

dydx=yy*0

r=sqrt(x^2+y^2)

ax=-v0^2*(x-r/re*(3./2*x*(x^2-y^2)/r^2+y*(2*x*y/r^2)))/(rc^2+x^2+y^2/q^2-r^3/re*(x^2-y^2)/r^2)

ay=-v0^2*(y/q^2-r/re*(3./2*y*(x^2+y^2)/r^2-x*(2*x*y/r^2)))/(rc^2+x^2+y^2/q^2-r^3/re*(x^2-y^2)/r^2)

dydx[0]=vx
dydx[1]=vy
dydx[2]=ax
dydx[3]=ay

return,dydx

end
