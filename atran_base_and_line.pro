; RK 9.5.2017
; return a baseline plus gaussian emission line based on an atran model to be fitted with mpfitfun
;RK 8.9.2017
; return the line flux corrected for the fitted transmission
;RK 12.10.2017 use Bill's algorithm,
;i.e. (baseline+line)*transmission convolved with instrument profile

function atran_base_and_line,wave,params
 
;  wave is the wavelength in microns, fixed step size monotinically increasing
;  parameters are:
;  p_alt    - the altidute in ft
;  p_pwv    -  the precipitable waver vapor in micron
;  p_za     - the zenith distance in degree 
;  p_cont   - continuum factor,
;  p_z      - zero level
;  p_flux   - line flux
;  p_center - line center
;  p_width  - line width
;  p_R      - Resolution (if non-positive or only one wavelength, then no smoothing)
  
  p_alt   =params[0]/1000.
  p_pwv   =params[1]
  p_za    =params[2]
  p_cont  =params[3]
  p_z     =params[4]
  p_flux  =params[5]
  p_center=params[6]
  p_width =params[7]
  p_R     =params[8]
  
  wave_and_trans=atran_transmission(wave[0],wave[-1],p_alt,p_pwv,p_za) ; full transmission spectrum
  wave_grid=wave_and_trans[*,0]
  trans=wave_and_trans[*,1]

  line =  fltarr(n_elements(wave_grid))                   ; creates the array for the emission line
  z=(wave_grid-p_center[0])/p_width[0]        ; argument for the gauss
  calc=where(abs(z) lt 4.29,cnt)             ; to improve cacluation time and avoid underflows (calculates the gaussian down to 1e-4)
  A0=p_flux[0]/(p_width[0]*sqrt(2.*!pi)) ; scaling of the flux argument
  if cnt gt 0 then line[calc]=A0*exp(-0.5*z[calc]^2)      ; Gaussian line

;flux correction if the width is too small for the sampling
  if n_elements(wave_grid) gt 1 then begin
     dlam=(wave_grid[-1]-wave_grid[0])/(n_elements(wave_grid)-1.)
     dummy=max(line,peak_idx)
     line[peak_idx] += p_flux[0]/dlam-total(line)
  endif

  
; multiply the continuum and line with the atmospheric transmisson
  spectrum=(p_cont+ line)*trans

  ;calculate the width of the kernel to smooth the model spectrum
  if (p_R gt 0 and n_elements(wave) gt 1) then begin
     lam=(wave_grid[-1]+wave_grid[0])/2.
     sigma=(lam/dlam)/p_R/2.35482

                                ;smooth the model spectrum and return the result
     spectrum=p_z+mygauss_smooth(spectrum,sigma,/edge_tru,/nan,$
                             width=min([6.*sigma,0.9*n_elements(spectrum)]))
  endif else begin
     spectrum+=p_z
  endelse

  return,interpol(spectrum,wave_grid,wave)
  
end
  
  
