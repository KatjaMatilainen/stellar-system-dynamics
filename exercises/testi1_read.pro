;***********************
; testi1_read.pro
; lukee tiedoston testi1
;*********************** 

openr,1,'testi1'
readf,1,a,b,c
readf,1,d,e,f
close,1

print,a,b,c,d,e,f
end