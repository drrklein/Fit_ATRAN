; RK 5.1.2017
; return an emission sprectrum based on an atran model to be fitted with mpfitfun

function atran_emission,wave,params
;  wave is the wavelength in microns
;  parameters are:
;  p_alt  - the altidute in ft
;  p_pwv  -  the precipitable waver vapor in micron
;  p_za   - the zenith distance in degree 
;  p_coef - a factor

  transmission=atran_baseline(wave,[params[0:2],1.,0.])

  return,params[3]*(1.-transmission)
end
  
  
