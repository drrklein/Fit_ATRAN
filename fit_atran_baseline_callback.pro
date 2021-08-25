;RK 4.1.2017
pro fit_atran_baseline_callback, status, error, bridge, userdata
  if status eq 2 then begin
     pwv_fit = bridge->getvar('pwv_fit')
     cont_fit = bridge->getvar('cont_fit')
     zero_fit = bridge->getvar('zero_fit')
     slope_fit = bridge->getvar('slope_fit')
     bl_subtracted = bridge->getvar('bl_subtracted')
     *(userdata.pout) = {status:status,pwv_fit:pwv_fit,cont_fit:cont_fit,zero_fit:zero_fit,slope_fit:slope_fit,$
                         bl_subtracted:bl_subtracted}
  endif else begin
     *(userdata.pout) = {status:status,instance:userdata.instance,error:error}
  endelse
end
