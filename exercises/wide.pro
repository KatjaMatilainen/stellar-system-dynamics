pro wide,window,exc=exc

;Hsalo ca 1995

if(keyword_set(exc)) then mexc=exc

if(n_params() eq 0) then begin

i=!d.window
newwin:


if(i ne -1) then begin

    if(keyword_set(exc)) then begin
        for j=0,n_elements(exc)-1 do begin
            if(i eq exc(j) and mexc(j) ne -1) then begin
                i=i-1
                mexc(j)=-1
                goto,skippi
            endif
        endfor
    endif

    wdelete,i
    i=!d.window
skippi:
    if(i gt 0) then goto, newwin
endif

endif else begin
    for j=0,n_elements(window)-1 do begin
        wdelete,window(j)
    endfor
endelse

end
