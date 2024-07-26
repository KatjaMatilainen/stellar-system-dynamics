;Aliohjelma, joka laskee rotaationopeuden eksponentiaaliselle kiekolle,
;kun sille syötetään vektori  R ja vakiot M ja Re

pro expdisc,R,M=M,Re=Re,vcirc=vcirc

;Käyttöohjeet

if (n_params() eq 0) then begin
   print,'pro expdisc,R,M,Re,vcirc'
   print,'R on vektori, M ja Re ovat vakioita'
endif

;Vakiot
G=1.d0
sig0=1.d0
Re=Re
R=R
M=M
const=4.*!dpi*G*sig0*Re

;y-vektori:
y=R/(2.*Re)
index=n_elements(R)
vcirc=findgen(index)*0.d0

;Analyyttinen kaava rotaationopeudelle:
for i=0,index-1 do begin
  vcirc(i)=sqrt(const*y(i)^2*(beseli(y(i),0)*beselk(y(i),0)-beseli(y(i),1)*$
  beselk(y(i),1)))
endfor

end
