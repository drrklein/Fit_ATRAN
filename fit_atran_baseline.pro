; RK 10.10.2016
; fit a baseline with an atran model
; RK 5.1.2017
; adding fitting an emission
; RK 12.5.2017
; changing the limits for zero and meaning of clip
; RK 15.5.2017
; adding a continuum slope and removing emission

pro fit_atran_baseline,cube,err_cube,header,base_idx,$
                       pwv_fit=pwv_fit,cont_fit=cont_fit,zero_fit=zero_fit,slope_fit=slope_fit,$
                       bl_subtracted=bl_subtracted,baselines=baselines,za=za,alt=alt,$
                       clip=clip,no_zero=no_zero,fix_pwv=fix_pwv,fix_slope=fix_slope,$
                       status=status,_extra=extra

s=size(cube)
wave=wave_hdr(header)

get_za_alt,header,za,alt

pwv_fit=fltarr(s[1],s[2])
cont_fit=fltarr(s[1],s[2])
zero_fit=fltarr(s[1],s[2])
slope_fit=fltarr(s[1],s[2])
bl_subtracted=cube

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0], relstep:0.1},7)

parinfo[0].fixed=1              ;altitude
if keyword_set(fix_pwv) then begin
   parinfo[1].fixed=1           ; pwv
   init_pwv=fix_pwv
endif else begin
   parinfo[1].limited=[1,1]     ; pwv
   parinfo[1].limits=[1,50]
   init_pwv=10
endelse
parinfo[2].fixed=1                ;za
parinfo[3].limited=[1,0]          ; cont
parinfo[3].limits=[0,0]         
if keyword_set(no_zero) then begin
   parinfo[4].fixed=1           ;zero
endif else begin
   if keyword_set(clip) then begin
      parinfo[4].limited=[1,1]  ; zero
      parinfo[4].limits=[-clip,+clip]
   endif
endelse

case n_elements(fix_slope) of   ; slope
   0:fix_slope=0.
   1:begin
      parinfo[5].fixed=1
   end
   2:begin
      parinfo[5].limited=1
      parinfo[5].limits=fix_slope
      fix_slope=total(fix_slope)/2.
   end
endcase

parinfo[6].fixed=1; wave0

for xp=0,s[1]-1 do begin
   for yp=0,s[2]-1 do begin
      good_base_idx=where(finite(cube[xp,yp,base_idx]))
      good_base_idx=base_idx[good_base_idx]
      if n_elements(good_base_idx) gt 10 then begin
         parinfo[*].value = [alt,init_pwv,za,median(cube[xp,yp,good_base_idx])>0,0,$
                             fix_slope,median(wave)]

         params=mpfitfun('atran_baseline',wave[good_base_idx],$
                         cube[xp,yp,good_base_idx],$
                         err_cube[xp,yp,good_base_idx],$
                         parinfo=parinfo,status=status,_extra=extra)
         pwv_fit[xp,yp]=params[1]
         cont_fit[xp,yp]=params[3]
         zero_fit[xp,yp]=params[4]
         slope_fit[xp,yp]=params[5]
         bl_subtracted[xp,yp,*]-=atran_baseline(wave,params)
      endif else begin
         pwv_fit[xp,yp]=!values.F_nan
         cont_fit[xp,yp]=!values.F_nan
         zero_fit[xp,yp]=!values.F_nan
         slope_fit[xp,yp]=!values.F_nan
         bl_subtracted[xp,yp,*]=!values.F_nan
      endelse
   endfor
   print,float(xp)/(s[1]-1)*100.,"% done"
endfor

baselines=cube-bl_subtracted


end
