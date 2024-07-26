pro plot_stamp,filename,header10=header10,scale=scale
if(n_params() le 0) then begin
    print,'plot_stamp,filename,header10=header10,scale=scale'
    print,'prints date+time and text given via filename'
    print,'HSalo 2003'
    return
endif




hostname=hostname()
whoami=mywhoami()
date=systime()
cd,curr=cwd
space='  '

header1=cwd+'/'+filename
header2=whoami+'@'+hostname+space+date
if(keyword_set(header10)) then header1=header10


sca=1.
if(keyword_set(scale)) then sca=scale
xyouts,0.01,0.005,header1,/nor,chars=.7*sca
xyouts,.99,0.005,header2,/nor,chars=.7*sca,ali=1
end
