;RK 12.10.2017
;get the transmission out of the atran hash (code used to be in atran_baseline.pro)

pro   atran_transmission_frac_idxs,values,p,val_idx,val_idx_f
  val_idx=(where(values-p le 0))[-1]
  
  if val_idx ne n_elements(values)-1 then begin
     val_idx_f=(float(p)-values[val_idx])/(values[val_idx+1]-values[val_idx])
  endif else begin
     val_idx-=1
     val_idx_f=1.
  endelse
end

function atran_transmission_old,wave,p_alt,p_pwv,p_za,smoothed=smoothed,interpolate=interpolate
;Not implemented, yet: if interpolate is not set, the nearest atran curve is used
  

; in the common block (smoothed and unsmoothed):
;  atran -hash with transmission: key is the orignial filename
;  wave_grid - wavelengths in micron for the transmissions
;  fn  - array with the filenames
;  alt - array with the sampled altitudes in 1000ft
;  pwv - array with the sampled pwv in micron
;  za  - array with the sampled zenith distances in deg
; smoothed_in_use - flag which set is in use
  
  common atran_transmission,atran_smoothed,wave_grid_smoothed,fn_smoothed,$
     alt_smoothed,pwv_smoothed,za_smoothed,$
     atran_unsmoothed,wave_grid_unsmoothed,fn_unsmoothed, wave_grid_base,$
     alt_unsmoothed,pwv_unsmoothed,za_unsmoothed,$
     atran,wave_grid,fn,alt,pwv,za,$
     smoothed_in_use

  if n_elements(smoothed_in_use) eq 0 then smoothed_in_use=-1 ; init
  
  if keyword_set(smoothed) then begin
     if n_elements(fn_smoothed) eq 0 then begin
        restore,'/Users/rklein1/IDL/Fit_ATRAN/atran_fifils_regrid.sav'
        atran_smoothed=atran
        wave_grid_smoothed=wave_grid
        fn_smoothed=fn
        alt_smoothed=alt
        pwv_smoothed=pwv
        za_smoothed=za
     endif
     if smoothed_in_use ne 1 then begin
        atran=atran_smoothed
        wave_grid=wave_grid_smoothed
        fn=fn_smoothed
        alt=alt_smoothed
        pwv=pwv_smoothed
        za=za_smoothed
        smoothed_in_use=1
     endif
  endif else begin
     if n_elements(fn_unsmoothed) eq 0 then begin
        restore,'/Users/rklein1/IDL/Fit_ATRAN/atran_fifils_regrid_unsmoothed.sav'
        atran_unsmoothed=atran
        wave_grid_unsmoothed=wave_grid
        fn_unsmoothed=fn
        alt_unsmoothed=alt
        pwv_unsmoothed=pwv
        za_unsmoothed=za
        wave_grid_base=(double(wave_grid[-1])/double(wave_grid[0]))^(1d0/(n_elements(wave_grid)-1)) ;log grid
     endif
     if smoothed_in_use ne 0 then begin
        atran=atran_unsmoothed
        wave_grid=wave_grid_unsmoothed
        fn=fn_unsmoothed
        alt=alt_unsmoothed
        pwv=pwv_unsmoothed
        za=za_unsmoothed
        smoothed_in_use=0
     endif        
  endelse

  atran_transmission_frac_idxs,alt,p_alt,alt_idx,alt_idx_f
  atran_transmission_frac_idxs,pwv,p_pwv,pwv_idx,pwv_idx_f
  atran_transmission_frac_idxs,za,p_za,za_idx,za_idx_f
  
  if keyword_set(smoothed) then begin
     wave_idx=interpol(findgen(n_elements(wave_grid)),wave_grid,wave)
  endif else begin
     wave_idx=alog(wave/double(wave_grid[0]))/alog(wave_grid_base)
  endelse

  key00=string(alt[alt_idx],ZA[za_idx],  pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key01=string(alt[alt_idx],ZA[za_idx],  pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key10=string(alt[alt_idx],ZA[za_idx+1],pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key11=string(alt[alt_idx],ZA[za_idx+1],pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')

  trans_4=[[[atran[key00]],[atran[key10]]],[[atran[key01]],[atran[key11]]]]
  trans=interpolate(trans_4,wave_idx,[za_idx_f],[pwv_idx_f],/grid)   

  if alt_idx_f ne 0 then begin
     key00=string(alt[alt_idx+1],ZA[za_idx],  pwv[pwv_idx]  ,$
                   format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
     key01=string(alt[alt_idx+1],ZA[za_idx],  pwv[pwv_idx+1],$
                   format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
     key10=string(alt[alt_idx+1],ZA[za_idx+1],pwv[pwv_idx]  ,$
                   format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
     key11=string(alt[alt_idx+1],ZA[za_idx+1],pwv[pwv_idx+1],$
                  format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')

     trans_4=[[[atran[key00]],[atran[key10]]],[[atran[key01]],[atran[key11]]]]
     trans*= (1-alt_idx_f)
     trans+=alt_idx_f*interpolate(trans_4,wave_idx,[za_idx_f],[pwv_idx_f],/grid)   
  endif
  return,trans
  end
  
  
