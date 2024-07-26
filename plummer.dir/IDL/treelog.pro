;treelog.pro
;convert TREELOG into TREELOGS
pro treelog,print=print,dire=dire


dir='./'
if(keyword_set(dire)) then dir=dire+'/'
print,'from: ',dir+'TREELOG   to:',dir+'TREELOGS'
close,1
openr,1,dir+'TREELOG'
a=''
readf,1,a
close,1

nlines=strlen(a)/81

nstep=-1
openw,1,dir+'TREELOGS'
for i=0,nlines-1 do begin
   b=strmid(a,i*81,81)
   if(keyword_set(print)) then print,b
   if(strpos(b,'time:') ne -1) then begin
      nstep=nstep+1
      if(nstep ge 0 and nstep lt 10)     then print,'SNAP00' +string(nstep,'(i1)')
      if(nstep ge 10 and nstep lt 100)   then print,'SNAP0'  +string(nstep,'(i2)')
      if(nstep ge 100 and nstep lt 1000) then print,'SNAP'   +string(nstep,'(i3)')
      print,b
   endif
   if(strpos(b,'ek') ne -1) then print,b
   if(strpos(b,'amx') ne -1) then print,b
   if(strpos(b,'cmpos') ne -1) then print,b
   
   printf,1,b
endfor
close,1
end
