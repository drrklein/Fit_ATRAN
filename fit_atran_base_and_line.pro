; RK 10.10.2016
; fit a baseline with an atran model
; RK 5.9.2017
; modifing to fil a line plus the baseline

pro decode_keyword_for_parinfo,xy_s,parinfo,parinfo_idx,keyword_val
    case n_elements(keyword_val) of
     0:case parinfo_idx of
        0: message,'alt should always be set'
        1: message,'pwv should always be set'
        2: message,'za should always be set'
        3: message,'cont should always be set'
        4: message,'zero should always be set'
        5: message,'flux should always be set'
        6: message,'need to set center'
        7: message,'need to set width'
        8: message,'resolution should always be set'
        else:message,'called with unknown parinfo_idx'
     endcase
     1:begin
        parinfo[parinfo_idx].fixed=1
        parinfo[parinfo_idx].value=keyword_val
     end        
     2:begin
        finite_kv=finite(keyword_val)
        finite_idx=where(finite_kv,cnt)
        parinfo[parinfo_idx].limited=finite_kv
        if cnt ne 0 then begin ;limited at least on one side or both
           parinfo[parinfo_idx].limits[finite_idx]=keyword_val[finite_idx]
           parinfo[parinfo_idx].value=total(keyword_val,/nan)/total(finite_kv)
        endif else begin ; unlimited
           parinfo[parinfo_idx].limits=[0.,0.]
           parinfo[parinfo_idx].value=0.
        endelse
        
     end
     3:begin
        finite_kv=finite(keyword_val)
        finite_idx=where(finite_kv[[0,2]],cnt)
        parinfo[parinfo_idx].limited=finite_kv[[0,2]]
        if cnt ne 0 then parinfo[parinfo_idx].limits[finite_idx]=keyword_val[([0,2])[finite_idx]]
        parinfo[parinfo_idx].value=keyword_val[1]
     end
     xy_s[2]:begin
        parinfo[parinfo_idx].fixed=1
        parinfo[parinfo_idx].value=keyword_val[xy_s[0],xy_s[1]]
     end
     else: message,'Elements of keyword value is unknown'
  endcase  
end

pro fit_atran_base_and_line,cube,err_cube,header,fit_wavel_idx,$
                            alt_fit=alt_fit,pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,$
                            flux_fit=flux_fit,center_fit=center_fit,width_fit=width_fit,res_fit=res_fit,$
                            alt_sig=alt_sig,pwv_sig=pwv_sig,cont_sig=cont_sig,zero_sig=zero_sig,$
                            flux_sig=flux_sig,center_sig=center_sig,width_sig=width_sig,res_sig=res_sig,$
                            residuals=residuals,base_and_lines=base_and_lines,za=za,fix_alt=fix_alt,$
                            fix_pwv=fix_pwv,fix_zero=fix_zero,fix_center=fix_center,$
                            fix_width=fix_width,fix_resolution=fix_resolution,status=status,$
                            fix_cont=fix_cont,fix_flux=fix_flux,$
                            _extra=extra

  s=size(cube)
  wave=wave_hdr(header)
  print,s
  get_za_alt,header,za,alt_temp
  if n_elements(fix_alt) eq 0 then begin
     fix_alt=alt_temp
  endif 
  if n_elements(fix_pwv) eq 0 then fix_pwv=[1.,7,50.]
  if n_elements(fix_resolution) eq 0 then fix_resolution=[500,1000.,2000]
  if n_elements(fix_cont) eq 0 then fix_cont=[0,median(cube)>0,!values.f_nan]
  if n_elements(fix_flux) eq 0 then fix_flux=[0,!values.f_nan]
  if n_elements(fix_zero) eq 0 then fix_zero=0
  
  alt_fit=fltarr(s[1],s[2])
  pwv_fit=fltarr(s[1],s[2])
  cont_fit=fltarr(s[1],s[2])
  zero_fit=fltarr(s[1],s[2])
  flux_fit=fltarr(s[1],s[2])
  center_fit=fltarr(s[1],s[2])
  width_fit=fltarr(s[1],s[2])
  res_fit=fltarr(s[1],s[2])
  alt_sig=fltarr(s[1],s[2])
  pwv_sig=fltarr(s[1],s[2])
  cont_sig=fltarr(s[1],s[2])
  zero_sig=fltarr(s[1],s[2])
  flux_sig=fltarr(s[1],s[2])
  center_sig=fltarr(s[1],s[2])
  width_sig=fltarr(s[1],s[2])
  res_sig=fltarr(s[1],s[2])
  residuals=cube
  
  
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0], step:0.1}, 9)
  
  
  for xp=0,s[1]-1 do begin
     for yp=0,s[2]-1 do begin

        xy_s=[xp,yp,s[1]*s[2]]
        decode_keyword_for_parinfo,xy_s,parinfo,0,fix_alt
        decode_keyword_for_parinfo,xy_s,parinfo,1,fix_pwv
        decode_keyword_for_parinfo,xy_s,parinfo,2,za
        decode_keyword_for_parinfo,xy_s,parinfo,3,fix_cont
        decode_keyword_for_parinfo,xy_s,parinfo,4,fix_zero
        decode_keyword_for_parinfo,xy_s,parinfo,5,fix_flux
        decode_keyword_for_parinfo,xy_s,parinfo,6,fix_center
        decode_keyword_for_parinfo,xy_s,parinfo,7,fix_width
        decode_keyword_for_parinfo,xy_s,parinfo,8,fix_resolution
        parinfo.step=[1000.,1.,5.,.1,.1,.1,$
                      parinfo[6].value/parinfo[8].value/3.,$ ;center
                      parinfo[6].value/parinfo[8].value/3.,$ ;width
                      100]
  
        good_idx=where(finite(cube[xp,yp,fit_wavel_idx]) and $
                       finite(err_cube[xp,yp,fit_wavel_idx]))
        wave_grid=wave[fit_wavel_idx[good_idx]]
        dlam=(wave_grid[-1]-wave_grid[0])/(n_elements(wave_grid)-1.)
        lam=(wave_grid[-1]+wave_grid[0])/2.
        sigma=(lam/dlam)/parinfo[8].value/2.35482

        if n_elements(wave_grid) ge $
           max([6.*sigma,n_elements(parinfo)-total(parinfo.fixed)]) $
        then begin
           
        
           params=mpfitfun('atran_base_and_line',wave_grid,$
                           cube[xp,yp,fit_wavel_idx[good_idx]],$
                           err_cube[xp,yp,fit_wavel_idx[good_idx]],$
                           parinfo=parinfo,status=status,errmsg=errmsg,$
                           _extra=extra,covar=covar,/quiet)
           if status le 0 then message,string(status)+errmsg
           
           alt_fit[xp,yp]=params[0]
           alt_sig[xp,yp]=sqrt(covar[0,0])
           pwv_fit[xp,yp]=params[1]
           pwv_sig[xp,yp]=sqrt(covar[1,1])
           cont_fit[xp,yp]=params[3]
           cont_sig[xp,yp]=sqrt(covar[3,3])
           zero_fit[xp,yp]=params[4]
           zero_sig[xp,yp]=sqrt(covar[4,4])
           flux_fit[xp,yp]=params[5]
           flux_sig[xp,yp]=sqrt(covar[5,5])
           center_fit[xp,yp]=params[6]
           center_sig[xp,yp]=sqrt(covar[6,6])
           width_fit[xp,yp]=params[7]
           width_sig[xp,yp]=sqrt(covar[7,7])
           res_fit[xp,yp]=params[8]
           res_sig[xp,yp]=sqrt(covar[8,8])
           residuals[xp,yp,*]-=atran_base_and_line(wave,params)
           
        endif else begin
           alt_fit[xp,yp]=!values.F_nan
           pwv_fit[xp,yp]=!values.F_nan
           cont_fit[xp,yp]=!values.F_nan
           zero_fit[xp,yp]=!values.F_nan
           flux_fit[xp,yp]=!values.F_nan
           center_fit[xp,yp]=!values.F_nan
           width_fit[xp,yp]=!values.F_nan
           res_fit[xp,yp]=!values.F_nan
           alt_sig[xp,yp]=!values.F_nan
           pwv_sig[xp,yp]=!values.F_nan
           cont_sig[xp,yp]=!values.F_nan
           zero_sig[xp,yp]=!values.F_nan
           flux_sig[xp,yp]=!values.F_nan
           center_sig[xp,yp]=!values.F_nan
           width_sig[xp,yp]=!values.F_nan
           res_sig[xp,yp]=!values.F_nan
           residuals[xp,yp,*]=!values.F_nan
        endelse
     endfor
     print,float(xp)/(s[1]-1)*100.,"% done"
  endfor
  
  base_and_lines=cube-residuals
  
  
end
