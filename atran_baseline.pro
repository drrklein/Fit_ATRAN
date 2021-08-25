; RK 9.10.2016
; return a baseline based on an atran model to be fitted with mpfitfun
; RK 15.5.2017
; adding a continuum slope

function atran_baseline,wave,params
  
;  wave is the wavelength in microns
;  parameters are:
;  p_alt  - the altidute in ft
;  p_pwv  -  the precipitable waver vapor in micron
;  p_za   - the zenith distance in degree 
;  p_cont - continuum factor,
;  p_z    - zero level
;  p_slope- continuum slope [per unit of wave, which is micron]
;  wave0  - wavelength around to pivot the continuum

  p_alt  =params[0]/1000.
  p_pwv  =params[1]
  p_za   =params[2]
  p_cont =params[3]
  p_z    =params[4]
  p_slope=params[5]
  wave0  =params[6]
  
  wave_and_trans = atran_transmission(wave[0],wave[-1],p_alt,p_pwv,p_za,/smoothed) 
  wave_grid=wave_and_trans[*,0]
  trans=wave_and_trans[*,1]

  transmission = interpol(trans,wave_grid,wave)
  
  return,(p_cont+(wave-wave0)*p_slope) * transmission + p_z
  
end
  
  
