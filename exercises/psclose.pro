;***************************************
;  PSCLOSE.PRO
;  closes ps-file
;  directs output to X

   pro psclose,junk=junk

;HS ca 1990
;***************************************

if(!d.name ne 'PS') then goto,fin

device,/CLOSE_FILE              ;,encapsul=0
!P.FONT=-1
if(keyword_set(junk)) then begin
    spawn,'lpr junk'
    print,'junk sent to printer'
endif

set_plot,'X'
fin:
end



