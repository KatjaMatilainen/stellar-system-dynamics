;***************************
;idlesim2.pro
;esimerkki IDl-aliohjelmasta
;ATK-kurssi 1992
;***************************

PRO ympyra,r,ala,piiri,print=print

if(n_params() eq 0) then begin
    print,'PRO ympyra,r,ala,piiri,print=print'
    return
endif

ala=!pi*r^2
piiri=2*!pi*r

if(keyword_set(print)) then begin
    print,'sade,piiri,ala'
    print,r,piiri,ala
endif

return
end
;***************************
