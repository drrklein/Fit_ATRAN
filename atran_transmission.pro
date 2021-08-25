;RK 12.10.2017
;get the transmission out of the atran hash (code used to be in
;atran_baseline.pro)

;RK 21.9.2018
;The returned transmission curve depended highly on the wavelength
;sampling, imnproved it by averaging the transmission curves to the
;sampling interval and then interpolating it

pro   atran_transmission_frac_idxs,values,p,val_idx,val_idx_f
  val_idx=(where(values-p le 0))[-1]
  
  if val_idx ne n_elements(values)-1 then begin
     val_idx_f=(float(p)-values[val_idx])/(values[val_idx+1]-values[val_idx])
  endif else begin
     val_idx-=1
     val_idx_f=1.
  endelse
end

function atran_transmission,wave_min,wave_max,p_alt,p_pwv,p_za,$
                            smoothed=smoothed



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
  
  key000=string(alt[alt_idx],ZA[za_idx],  pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key001=string(alt[alt_idx],ZA[za_idx],  pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key010=string(alt[alt_idx],ZA[za_idx+1],pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key011=string(alt[alt_idx],ZA[za_idx+1],pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key100=string(alt[alt_idx+1],ZA[za_idx],  pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key101=string(alt[alt_idx+1],ZA[za_idx],  pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key110=string(alt[alt_idx+1],ZA[za_idx+1],pwv[pwv_idx]  ,$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
  key111=string(alt[alt_idx+1],ZA[za_idx+1],pwv[pwv_idx+1],$
               format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')

  trans_8=[[[atran[key000]],[atran[key010]]],[[atran[key001]],[atran[key011]]],$
           [[atran[key100]],[atran[key110]]],[[atran[key101]],[atran[key111]]]]

  trans_8=reform(trans_8,n_elements(wave_grid),2,2,2,/overwrite)
  trans=interpolate(trans_8,[alt_idx_f],[za_idx_f],[pwv_idx_f],/grid)   
 
  idx=where(wave_min le wave_grid and wave_grid le wave_max)
  idx=[idx[0]-1,idx,idx[-1]+1]
  
  return,[[wave_grid[idx]],[trans[idx]]]
  end
  
  
