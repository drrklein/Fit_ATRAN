;RK 4.1.2017
;modified for a general applicion
pro callback, status, error, bridge, userdata
  if status eq 2 then begin
     out = bridge->getvar('out')
     *(userdata.pout) = out
  endif else begin
     *(userdata.pout) = {status:status,error:error}
  endelse
end
