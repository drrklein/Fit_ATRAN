;RK 9.5.2017
pro fit_atran_base_and_line_callback, status, error, bridge, userdata
  if status eq 2 then begin
     alt_fit = bridge->getvar('alt_fit')
     pwv_fit = bridge->getvar('pwv_fit')
     cont_fit = bridge->getvar('cont_fit')
     zero_fit = bridge->getvar('zero_fit')
     flux_fit = bridge->getvar('flux_fit')
     center_fit = bridge->getvar('center_fit')
     width_fit = bridge->getvar('width_fit')
     res_fit = bridge->getvar('res_fit')
     alt_sig = bridge->getvar('alt_sig')
     pwv_sig = bridge->getvar('pwv_sig')
     cont_sig = bridge->getvar('cont_sig')
     zero_sig = bridge->getvar('zero_sig')
     flux_sig = bridge->getvar('flux_sig')
     center_sig = bridge->getvar('center_sig')
     width_sig = bridge->getvar('width_sig')
     res_sig = bridge->getvar('res_sig')
     residuals = bridge->getvar('residuals')
     
     *(userdata.pout) = {status:status,alt_fit:alt_fit,pwv_fit:pwv_fit,cont_fit:cont_fit,$
                         zero_fit:zero_fit,flux_fit:flux_fit,center_fit:center_fit,$
                         width_fit:width_fit,res_fit:res_fit,alt_sig:alt_sig,$
                         pwv_sig:pwv_sig,cont_sig:cont_sig,$
                         zero_sig:zero_sig,flux_sig:flux_sig,center_sig:center_sig,$
                         width_sig:width_sig,res_sig:res_sig,residuals:residuals}
  endif else begin
     *(userdata.pout) = {status:status,instance:userdata.instance,error:error}
  endelse
end
